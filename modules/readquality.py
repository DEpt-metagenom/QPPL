import os
import subprocess
import re
import logging
import time
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

# Function to generate commands for read quality tasks
def generate_commands(read_quality_dir, input_dir, file, basename, task_params):
    return [
        f"echo \"Tool version (seqkit):\" $(seqkit version)",
        f"seqkit stat -a -j {task_params['readquality']['seqkit_threads']} -T {input_dir}/{file} > {read_quality_dir}/{basename}_seqkit.tsv",
        f"echo \"Tool version (NanoQC):\" $(nanoQC --version)",
        f"nanoQC -o {read_quality_dir}/ {input_dir}/{file}",
        f"echo \"Tool version (nanoq):\" $(nanoq --version)",
        f"nanoq -i {input_dir}/{file} -s -H > {read_quality_dir}/{basename}_nanoq.log"
    ]

# Function to move files if they exist
def move_file_if_exists(src, dst):
    if os.path.isfile(src):
        shutil.move(src, dst)
        logger.info(f"Moved {src} to {dst}")
    else:
        logger.warning(f"{src} not found.")

# Main function to run read quality tasks
def run_read_quality(input_dir, output_dir, input_files, task_params):
    read_quality_dir = os.path.join(output_dir, "readquality_task")

    # Ensure the read quality directory exists
    os.makedirs(read_quality_dir, exist_ok=True)
    log_path = os.path.join(output_dir, "readquality.log")
    setup_logging(log_path)

    logger.info(f"Using folder: \"{read_quality_dir}\"")

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
        logger.info(f"Starting read quality tasks for file: {file}")

        # Run commands in environments
        for step_number, env in envs.items():
            commands = generate_commands(read_quality_dir, input_dir, file, basename, task_params)
            for command in commands:
                run_command(env["env_name"], command, file, step_number)

        # Move log files if they exist
        move_file_if_exists(os.path.join(read_quality_dir, "NanoQC.log"), os.path.join(read_quality_dir, f"{basename}_NanoQC.log"))
        move_file_if_exists(os.path.join(read_quality_dir, "nanoQC.html"), os.path.join(read_quality_dir, f"{basename}_NanoQC.html"))

        # Calculate elapsed time per file
        elapsed = time.time() - start_time
        logger.info(f"Finished read quality tasks for file: {file} (elapsed time: {elapsed:.2f} seconds)")

    # Calculate total elapsed time
    total_elapsed = time.time() - total_start_time
    logger.info(f"Total elapsed time for all tasks: {total_elapsed:.2f} seconds")