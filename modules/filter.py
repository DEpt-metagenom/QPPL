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

    # Resolve which noadapter file already exists (uncompressed takes priority since
    # it is the freshest porechop output; .gz is the archived form from a prior run).
    noadapter_fastq = output_file_porechop
    noadapter_gz    = output_file_porechop + ".gz"
    noadapter_cached = os.path.exists(noadapter_fastq) or os.path.exists(noadapter_gz)

    # Add commands based on whether to cut adapters
    if task_params['filter']['cut_adapter'] == 'y':
        if noadapter_cached:
            logger.info(f"Skipping porechop for {file} — cached noadapter file found")
        else:
            commands.append(f"echo \"Tool version (porechop_abi):\" $(porechop_abi --version)")
            commands.append(f"porechop_abi -abi --discard_middle -i {input_dir}/{file} -o {output_file_porechop} --format fastq -t {task_params['filter']['porechop_threads']} --verbosity 0")
            commands.append(f"gzip {output_file_porechop}")
            commands.append(f"rm -r tmp")

        # Point to whichever noadapter form exists; prefer uncompressed (no gunzip overhead).
        if os.path.exists(noadapter_fastq):
            input_file = noadapter_fastq
            read_cmd = f"cat {input_file}"
        else:
            input_file = noadapter_gz
            read_cmd = f"gunzip -c {input_file}"
    else:
        input_file = os.path.join(input_dir, file)
        read_cmd = f"gunzip -c {input_file}"

    # Command for Seqkit, NanoLyse and NanoFilt
    commands.append(f"echo \"Tool version (Seqkit):\" $(seqkit version)")
    commands.append(f"seqkit stat -a -j {task_params['filter']['seqkit_threads']} -T {input_dir}/{file} > {filter_dir}/{basename}_pre_tmp.tsv")
    commands.append(
        f"""awk 'BEGIN{{FS=OFS="\\t"}} NR==1{{for(i=1;i<=NF;i++) if($i!="file" && $i!="format" && $i!="type") $i=$i"_PRE"}} 1' {filter_dir}/{basename}_pre_tmp.tsv > {filter_dir}/{basename}_seqkit_PRE.tsv"""
    )
    commands.append(f"echo \"Tool version (NanoLyse):\" $(NanoLyse --version)")
    commands.append(f"echo \"Tool version (NanoFilt):\" $(NanoFilt --version)")
    commands.append(f"{read_cmd} | NanoLyse -r {filter_dir}/lambda.fasta | NanoLyse -r {filter_dir}/DNA_CS.fasta | NanoFilt -q {task_params['filter']['filter_quality']} -l {task_params['filter']['filter_length']} | gzip > {output_file_nanotools}")
    commands.append(f"seqkit stat -a -j {task_params['filter']['seqkit_threads']} -T {output_file_nanotools} > {filter_dir}/{basename}_post_tmp.tsv")
    commands.append(
        f"""awk 'BEGIN{{FS=OFS="\\t"}} NR==1{{for(i=1;i<=NF;i++) if($i!="file" && $i!="format" && $i!="type") $i=$i"_POST"}} 1' {filter_dir}/{basename}_post_tmp.tsv > {filter_dir}/{basename}_seqkit_POST.tsv"""
    )
    commands.append(
        f"""awk 'BEGIN{{FS=OFS="\\t"}}
        function norm(x, y) {{
            y=x
            sub(/^.*\\//, "", y)
            sub(/^filtered_/, "", y)
            sub(/^noadapter_/, "", y)
            sub(/\\.(fastq|fq)(\\.gz)?$/, "", y)
            return y
        }}
        NR==FNR {{
            if (FNR==1) {{for(i=4;i<=NF;i++) pre_h[i]=$i; next}}
            key=norm($1) OFS $2 OFS $3
            for(i=4;i<=NF;i++) pre[key,i]=$i
            next
        }}
        FNR==1 {{
            header="file" OFS "format" OFS "type"
            for(i=4;i<=NF;i++) header=header OFS pre_h[i] OFS $i
            print header
            next
        }}
        {{
            key=norm($1) OFS $2 OFS $3
            line=norm($1) OFS $2 OFS $3
            for(i=4;i<=NF;i++) line=line OFS (((key SUBSEP i) in pre)?pre[key,i]:"NA") OFS $i
            print line
        }}' {filter_dir}/{basename}_seqkit_PRE.tsv {filter_dir}/{basename}_seqkit_POST.tsv > {filter_dir}/{basename}_seqkit_PRE_POST.tsv"""
    )
    commands.append(f"rm NanoLyse.log {filter_dir}/{basename}_pre_tmp.tsv {filter_dir}/{basename}_post_tmp.tsv {filter_dir}/{basename}_seqkit_PRE.tsv {filter_dir}/{basename}_seqkit_POST.tsv")

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
        done_marker = os.path.join(filter_dir, f".{basename}.done")

        if os.path.exists(done_marker):
            logger.info(f"Skipping file: {file} (already filtered in a previous run)")
            continue

        # Start task timer per file
        start_time = time.time()

        logger.info("=" * 60)
        logger.info(f"Starting filter tasks for file: {file}")

        # Run commands in environments
        for step_number in sorted(envs.keys()):
            commands = generate_commands(input_dir, filter_dir, file, task_params, basename)
            for command in commands:
                run_command(envs[step_number]["env_name"], command, file, step_number)

        open(done_marker, "w").close()

        # Calculate elapsed time per file
        elapsed = time.time() - start_time
        logger.info(f"Finished filter tasks for file: {file} (elapsed time: {elapsed:.2f} seconds)")

    # Calculate total elapsed time
    total_elapsed = time.time() - total_start_time
    logger.info(f"Total elapsed time for all tasks: {total_elapsed:.2f} seconds")