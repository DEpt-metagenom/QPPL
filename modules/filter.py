import os
import subprocess
import re
import logging
import time
import requests
import gzip
import shutil

# Call logger
logger = logging.getLogger("QPPL")

# Function to fully log messages only to file
def log_to_file_only(message, level=logging.INFO):
    for h in logger.handlers:
        if isinstance(h, logging.FileHandler):
            record = logger.makeRecord(logger.name, level, None, None, message, None, None)
            h.emit(record)

# Function to set up logging
def setup_logging(log_path):
    file_handler = logging.FileHandler(log_path)
    file_handler.setFormatter(logging.Formatter("[%(levelname)s] %(asctime)s - %(message)s", "%Y-%m-%d %H:%M:%S"))

    # Remove previous file handlers to avoid duplicate logs
    for h in logger.handlers[:]:
        if isinstance(h, logging.FileHandler):
            logger.removeHandler(h)
    logger.addHandler(file_handler)

# Function to run a command in a conda environment
def run_command(env_name, command, file, env_number):
    log_to_file_only(f"Running command in {env_name} on {file}: $ {command}", logging.INFO)
    result = subprocess.run(["conda", "run", "-n", env_name, "bash", "-c", command], text=True, capture_output=True)

    if result.returncode != 0:
        logger.error(f"Tool failed for file: {file} (step {env_number})\nCommand: {command}\nError: {result.stderr.strip()}")
        return
    
    out = result.stdout.strip()
    if not out:
        return

    is_version_line = "Tool version" in command
    log_fn = logger.debug if is_version_line else logger.info
    for line in out.splitlines():
        log_fn(line)

# Constants for URLs
LAMBDA_FASTA_URL = "https://github.com/wdecoster/nanolyse/raw/master/reference/lambda.fasta.gz"
DNA_CS_FASTA_URL = "https://raw.githubusercontent.com/wdecoster/nanolyse/master/reference/DNA_CS.fasta"

# Function to download a file from a URL
def download_file(url, dest):
    response = requests.get(url, stream=True)
    if response.status_code == 200:
        with open(dest, 'wb') as f:
            f.write(response.content)
    else:
        logger.error(f"Failed to download {url}: {response.status_code}")

# Function to gunzip a file
def gunzip_file(gz_file, dest):
    with gzip.open(gz_file, 'rb') as f_in:
        with open(dest, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

# Function to generate commands for filtering tasks
def generate_commands(input_dir, filter_dir, file, task_params, basename):
    output_file_porechop = os.path.join(filter_dir, f"noadapter_{basename}.fastq")
    output_file_nanotools = os.path.join(filter_dir, f"filtered_{file}")

    commands = []

    # Download reference files if they don't exist
    lambda_fasta_path = os.path.join(filter_dir, "lambda.fasta")
    if not os.path.exists(lambda_fasta_path):
        lambda_gz_path = lambda_fasta_path + ".gz"
        download_file(LAMBDA_FASTA_URL, lambda_gz_path)
        gunzip_file(lambda_gz_path, lambda_fasta_path)

    dna_cs_fasta_path = os.path.join(filter_dir, "DNA_CS.fasta")
    if not os.path.exists(dna_cs_fasta_path):
        download_file(DNA_CS_FASTA_URL, dna_cs_fasta_path)

    # Add commands based on whether to cut adapters
    if task_params['filter']['cut_adapter'] == 'y':
        commands.append(f"echo \"Tool version (porechop_abi):\" $(porechop_abi --version)")
        commands.append(f"porechop_abi -abi --discard_middle -i {input_dir}/{file} -o {output_file_porechop} --format fastq -t {task_params['filter']['porechop_threads']} --verbosity 0")
        commands.append(f"gzip {output_file_porechop}")
        input_file = output_file_porechop
        commands.append(f"rm -r tmp")
    else:
        input_file = os.path.join(input_dir, file)

    # Command for NanoLyse and NanoFilt
    commands.append(f"echo \"Tool version (NanoLyse):\" $(NanoLyse --version)")
    commands.append(f"echo \"Tool version (NanoFilt):\" $(NanoFilt --version)")
    commands.append(f"gunzip -c {input_file} | NanoLyse -r {filter_dir}/lambda.fasta | NanoLyse -r {filter_dir}/DNA_CS.fasta | NanoFilt -q {task_params['filter']['filter_quality']} -l {task_params['filter']['filter_length']} | gzip > {output_file_nanotools}")
    commands.append(f"rm NanoLyse.log")

    return commands

# Main function to run filter tasks
def run_filter(input_dir, output_dir, input_files, task_params):
    filter_dir = os.path.join(output_dir, "filter_task")

    # Ensure the filter directory exists
    os.makedirs(filter_dir, exist_ok=True)
    log_path = os.path.join(output_dir, "filter.log")
    setup_logging(log_path)

    logger.info(f"Using folder: \"{filter_dir}\"")

    # Start total timer
    total_start_time = time.time()

    # Define used environments
    envs = {
        1: {"env_name": "QPPL_env"}
    }

    # Process each input file
    for file in input_files:
        basename = re.sub(r'(\.fastq|\.fastq\.gz|\.fq|\.fq\.gz)$', '', file, flags=re.IGNORECASE)
        
        # Start task timer per file
        start_time = time.time()

        logger.info("=" * 60)
        logger.info(f"Starting filter tasks for file: {file}")

        # Run commands in environments
        for step_number in sorted(envs.keys()):
            commands = generate_commands(input_dir, filter_dir, file, task_params, basename)
            for command in commands:
                run_command(envs[step_number]["env_name"], command, file, step_number)

        # Calculate elapsed time per file
        elapsed = time.time() - start_time
        logger.info(f"Finished filter tasks for file: {file} (elapsed time: {elapsed:.2f} seconds)")

    # Calculate total elapsed time
    total_elapsed = time.time() - total_start_time
    logger.info(f"Total elapsed time for all tasks: {total_elapsed:.2f} seconds")