import os
import subprocess
import re
import logging
import time

# Call logger
logger = logging.getLogger("QPPL")

# Output subdirectories for assembly task inside assembly directory
OUTPUT_DIRS = {
    1: "out_flye",
    2: "out_goldrush",
    3: "out_raven",
    4: "out_shasta",
    5: "out_wtdbg2",
    6: "out_polish",
    7: "out_checkv",
    8: "contigs",
    9: "out_vclust",
    10: "depth_coverage"
}

# All assemblers used in the assembly task
ASSEMBLERS = {
    1: "flye",
    2: "goldrush",
    3: "raven",
    4: "shasta",
    5: "wtdbg2"
}

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

# Function to generate commands for assembly tasks per environment
def generate_commands(assembly_dir, output_dir, file, task_params, basename, step_number):
    meta_flag = "--meta" if task_params['assembly']['flye_meta'] == 'y' else ""
    wtdbg2_realign_flag = "-R" if task_params['assembly']['wtdbg2_realign'] == 'y' else ""

    contig_prefixes = {
        1: "flye_contig_{nr}",
        2: "goldrush_contig_{nr}",
        3: "raven_contig_{nr}",
        4: "shasta_contig_{nr}",
        5: "wtdbg2_contig_{nr}",
    }

    score_formula = "((1-(host_genes/gene_count))*(viral_genes/(contig_length*kmer_freq))*(completeness/100)*((100-contamination)/100)*1000)"

    def create_subsubdir(dir_path):
        return f"mkdir -p {dir_path}"

    def seqkit_replace_command(prefix, input_file, output_file):
        return f"seqkit replace -p '.+' -r {prefix} {input_file} | seqkit seq -w 0 > {output_file}"

    def remove_subsubdir(dir_path):
        return f"rm -r {dir_path}"

    commands_by_step = {
        0: [ # Gunzipping input files
            f"gunzip -c {output_dir}/filter_task/filtered_{file} > {assembly_dir}/filtered_{basename}.fastq"
        ],
        1: [  # Flye
            f"echo \"Tool version (Flye):\" $(flye --version)",
            f"echo \"Tool version with Flye (seqkit):\" $(seqkit version)",
            f"flye --{task_params['assembly']['flye_mode']} {assembly_dir}/filtered_{basename}.fastq -o {assembly_dir}/{OUTPUT_DIRS[1]}/{basename} --threads {task_params['assembly']['flye_threads']} {meta_flag}",
            seqkit_replace_command(contig_prefixes[1], f"{assembly_dir}/{OUTPUT_DIRS[1]}/{basename}/assembly.fasta", f"{assembly_dir}/{OUTPUT_DIRS[1]}/{basename}_flye.fasta"),
            remove_subsubdir(f"{assembly_dir}/{OUTPUT_DIRS[1]}/{basename}")
        ],
        2: [ # Goldrush
            f"echo \"Tool version (goldrush):\" $(goldrush)",
            f"echo \"Tool version with Goldrush (seqkit):\" $(seqkit version)",
            create_subsubdir(f"{assembly_dir}/{OUTPUT_DIRS[2]}/{basename}"),
            f"cp {assembly_dir}/filtered_{basename}.fastq {assembly_dir}/{OUTPUT_DIRS[2]}/{basename}/tmp.fq",
            f"cd {assembly_dir}/{OUTPUT_DIRS[2]}/{basename} && goldrush run reads=tmp p=tmp G={task_params['assembly']['goldrush_genome_size']} m={task_params['assembly']['goldrush_min_read_length']} t={task_params['assembly']['goldrush_threads']} && cd -",
            f"mv {assembly_dir}/{OUTPUT_DIRS[2]}/{basename}/goldrush_intermediate_files/tmp*.polished.fa {assembly_dir}/{OUTPUT_DIRS[2]}/{basename}/assembly.fasta",
            seqkit_replace_command(contig_prefixes[2], f"{assembly_dir}/{OUTPUT_DIRS[2]}/{basename}/assembly.fasta", f"{assembly_dir}/{OUTPUT_DIRS[2]}/{basename}_goldrush.fasta"),
            remove_subsubdir(f"{assembly_dir}/{OUTPUT_DIRS[2]}/{basename}")
        ],
        3: [ # Raven
            f"echo \"Tool version (raven):\" $(raven --version)",
            f"echo \"Tool version with Raven (seqkit):\" $(seqkit version)",
            create_subsubdir(f"{assembly_dir}/{OUTPUT_DIRS[3]}/{basename}"),
            f"raven --disable-checkpoints -t {task_params['assembly']['raven_threads']} {assembly_dir}/filtered_{basename}.fastq > {assembly_dir}/{OUTPUT_DIRS[3]}/{basename}/{basename}_temp.fasta",
            seqkit_replace_command(contig_prefixes[3], f"{assembly_dir}/{OUTPUT_DIRS[3]}/{basename}/{basename}_temp.fasta", f"{assembly_dir}/{OUTPUT_DIRS[3]}/{basename}_raven.fasta"),
            remove_subsubdir(f"{assembly_dir}/{OUTPUT_DIRS[3]}/{basename}"),
        ],
        4: [ # Shasta
            f"echo \"Tool version (shasta):\" $(shasta --version)",
            f"echo \"Tool version with Shasta (seqkit):\" $(seqkit version)",
            f"seqkit grep -v -s -r -p '[^ACGT]' {assembly_dir}/filtered_{basename}.fastq > {assembly_dir}/{OUTPUT_DIRS[4]}/{basename}_clean.fq",
            f"shasta --input {assembly_dir}/{OUTPUT_DIRS[4]}/{basename}_clean.fq --config Nanopore-May2022 --threads {task_params['assembly']['shasta_threads']} --Reads.minReadLength {task_params['assembly']['shasta_read_length']} --assemblyDirectory {assembly_dir}/{OUTPUT_DIRS[4]}/{basename}",
            seqkit_replace_command(contig_prefixes[4], f"{assembly_dir}/{OUTPUT_DIRS[4]}/{basename}/Assembly.fasta", f"{assembly_dir}/{OUTPUT_DIRS[4]}/{basename}_shasta.fasta"),
            f"rm {assembly_dir}/{OUTPUT_DIRS[4]}/{basename}_clean.fq",
            remove_subsubdir(f"{assembly_dir}/{OUTPUT_DIRS[4]}/{basename}")
        ],
        5: [ # WTDBG2
            f"echo \"Tool version (wtdbg2):\" $(wtdbg2 --version)",
            f"echo \"Tool version with WTDBG2 (seqkit):\" $(seqkit version)",
            create_subsubdir(f"{assembly_dir}/{OUTPUT_DIRS[5]}/{basename}"),
            f"wtdbg2 -g {task_params['assembly']['wtdbg2_genome_size']} -t {task_params['assembly']['wtdbg2_threads']} -i {assembly_dir}/filtered_{basename}.fastq -fo {assembly_dir}/{OUTPUT_DIRS[5]}/{basename}/{basename} -L {task_params['assembly']['wtdbg2_read_length']} -k {task_params['assembly']['wtdbg2_kmer_fsize']} -p {task_params['assembly']['wtdbg2_kmer_psize']} -s {task_params['assembly']['wtdbg2_min_similarity']} {wtdbg2_realign_flag}",
            f"wtpoa-cns -t {task_params['assembly']['wtdbg2_threads']} -i {assembly_dir}/{OUTPUT_DIRS[5]}/{basename}/{basename}.ctg.lay.gz -fo {assembly_dir}/{OUTPUT_DIRS[5]}/{basename}/{basename}.ctg.fasta",
            seqkit_replace_command(contig_prefixes[5], f"{assembly_dir}/{OUTPUT_DIRS[5]}/{basename}/{basename}.ctg.fasta", f"{assembly_dir}/{OUTPUT_DIRS[5]}/{basename}_wtdbg2.fasta"),
            remove_subsubdir(f"{assembly_dir}/{OUTPUT_DIRS[5]}/{basename}")
        ],
        6: [ # Polishing
            f"echo \"Tool version (minimap2):\" $(minimap2 --version)",
            f"echo \"Tool version (racon):\" $(racon --version)",
            f"echo \"Tool version (medaka):\" $(medaka --version)",
            f"cp {assembly_dir}/filtered_{basename}.fastq {assembly_dir}/{OUTPUT_DIRS[6]}/tmp_{basename}.fastq",
            *[rf"for asm in {' '.join(ASSEMBLERS.values())}; do minimap2 -x map-ont -t {task_params['assembly']['polish_threads']} -a {assembly_dir}/out_$asm/{basename}_$asm.fasta {assembly_dir}/{OUTPUT_DIRS[6]}/tmp_{basename}.fastq > {assembly_dir}/{OUTPUT_DIRS[6]}/{basename}_align_$asm.sam; done"],
            *[rf"for asm in {' '.join(ASSEMBLERS.values())}; do racon -t {task_params['assembly']['polish_threads']} -f {assembly_dir}/{OUTPUT_DIRS[6]}/tmp_{basename}.fastq {assembly_dir}/{OUTPUT_DIRS[6]}/{basename}_align_$asm.sam {assembly_dir}/out_$asm/{basename}_$asm.fasta > {assembly_dir}/{OUTPUT_DIRS[6]}/{basename}_asmbl_$asm.fa; done"],
            *[rf"for asm in {' '.join(ASSEMBLERS.values())}; do mini_align -i {assembly_dir}/{OUTPUT_DIRS[6]}/tmp_{basename}.fastq -r {assembly_dir}/{OUTPUT_DIRS[6]}/{basename}_asmbl_$asm.fa -P -m -p {assembly_dir}/{OUTPUT_DIRS[6]}/{basename}_align_$asm -t {task_params['assembly']['polish_threads']}; done"],
            *[rf"for asm in {' '.join(ASSEMBLERS.values())}; do medaka inference --model r1041_e82_400bps_sup_g615 {assembly_dir}/{OUTPUT_DIRS[6]}/{basename}_align_$asm.bam {assembly_dir}/{OUTPUT_DIRS[6]}/{basename}_align_$asm.hdf; done"],
            *[rf"for asm in {' '.join(ASSEMBLERS.values())}; do medaka sequence {assembly_dir}/{OUTPUT_DIRS[6]}/{basename}_align_$asm.hdf {assembly_dir}/{OUTPUT_DIRS[6]}/{basename}_asmbl_$asm.fa {assembly_dir}/{OUTPUT_DIRS[6]}/{basename}_polished_$asm.fa; done"],
            rf"find {assembly_dir}/{OUTPUT_DIRS[6]}/ -type f ! -name '*_polished_*.fa' -exec rm {{}} \;"
        ],
        7: [ # Quality assessment
            f"echo \"Tool version (CheckV):\" $(checkv -h 2>&1 | head -n1 | cut -d':' -f1)",
            rf"mkdir {assembly_dir}/{OUTPUT_DIRS[7]}/{basename}",
            *[rf"for asm in {' '.join(ASSEMBLERS.values())}; do cat {assembly_dir}/{OUTPUT_DIRS[6]}/{basename}_polished_$asm.fa >> {assembly_dir}/{OUTPUT_DIRS[7]}/{basename}/{basename}_polished.fa; done"],
            f"checkv end_to_end {assembly_dir}/{OUTPUT_DIRS[7]}/{basename}/{basename}_polished.fa {assembly_dir}/{OUTPUT_DIRS[7]}/{basename} -d {task_params['assembly']['checkv_db']} -t {task_params['assembly']['checkv_threads']}"
        ],
        8: [ # Low or medium quality and contaminated contig removal
            rf"csvcut -t -c contig_id,contig_length,gene_count,viral_genes,host_genes,checkv_quality,miuvig_quality,completeness,completeness_method,contamination,kmer_freq {assembly_dir}/{OUTPUT_DIRS[7]}/{basename}/quality_summary.tsv | csvgrep -c miuvig_quality -m High-quality | awk -F, '$5 <= $4' > {assembly_dir}/{OUTPUT_DIRS[7]}/{basename}/{basename}_summary.csv",
            rf"mkdir {assembly_dir}/{OUTPUT_DIRS[8]}/{basename}",
            rf"awk -F',' 'NR>1 {{print $1}}' {assembly_dir}/{OUTPUT_DIRS[7]}/{basename}/{basename}_summary.csv | while IFS= read -r contig; do seqtk subseq {assembly_dir}/{OUTPUT_DIRS[7]}/{basename}/{basename}_polished.fa <(echo $contig) > {assembly_dir}/{OUTPUT_DIRS[8]}/{basename}/$contig.fa; done"
        ],
        9: [ # Clustering and selection of representative contigs
            f"echo \"Tool version (vclust):\" $(vclust --version)",
            rf"mkdir {assembly_dir}/{OUTPUT_DIRS[9]}/{basename}",
            rf"cat {assembly_dir}/{OUTPUT_DIRS[8]}/{basename}/*.fa > {assembly_dir}/{OUTPUT_DIRS[9]}/{basename}/{basename}_contigs.fa",
            rf"if [ $(grep -c '>' {assembly_dir}/{OUTPUT_DIRS[9]}/{basename}/{basename}_contigs.fa) -eq 1 ]; then cp {assembly_dir}/{OUTPUT_DIRS[9]}/{basename}/{basename}_contigs.fa {assembly_dir}/{OUTPUT_DIRS[9]}/{basename}/{basename}_selected_genomes.fasta; fi",
            rf"if [ ! -f {assembly_dir}/{OUTPUT_DIRS[9]}/{basename}/{basename}_selected_genomes.fasta ]; then vclust prefilter -i {assembly_dir}/{OUTPUT_DIRS[9]}/{basename}/{basename}_contigs.fa -o {assembly_dir}/{OUTPUT_DIRS[9]}/{basename}/{basename}_fltr.txt --min-ident 0.7; fi",
            rf"if [ ! -f {assembly_dir}/{OUTPUT_DIRS[9]}/{basename}/{basename}_selected_genomes.fasta ]; then vclust align -i {assembly_dir}/{OUTPUT_DIRS[9]}/{basename}/{basename}_contigs.fa -o {assembly_dir}/{OUTPUT_DIRS[9]}/{basename}/{basename}_align.tsv --filter {assembly_dir}/{OUTPUT_DIRS[9]}/{basename}/{basename}_fltr.txt --outfmt lite --out-ani 0.95 --out-qcov 0.85; fi",
            rf"if [ ! -f {assembly_dir}/{OUTPUT_DIRS[9]}/{basename}/{basename}_selected_genomes.fasta ]; then vclust cluster -i {assembly_dir}/{OUTPUT_DIRS[9]}/{basename}/{basename}_align.tsv -o {assembly_dir}/{OUTPUT_DIRS[9]}/{basename}/{basename}_clusters.tsv --ids {assembly_dir}/{OUTPUT_DIRS[9]}/{basename}/{basename}_align.ids.tsv --algorithm complete --metric ani --ani 0.95 --qcov 0.85; fi",
            rf"if [ ! -f {assembly_dir}/{OUTPUT_DIRS[9]}/{basename}/{basename}_selected_genomes.fasta ]; then awk 'NR>1 {{print $1,$2}}' {assembly_dir}/{OUTPUT_DIRS[9]}/{basename}/{basename}_clusters.tsv | sort -k1,1 | tr ' ' ',' > {assembly_dir}/{OUTPUT_DIRS[9]}/{basename}/{basename}_sorted_clusters.csv; else awk 'NR==1 {{print substr($0,2),0}}' {assembly_dir}/{OUTPUT_DIRS[9]}/{basename}/{basename}_selected_genomes.fasta > {assembly_dir}/{OUTPUT_DIRS[9]}/{basename}/{basename}_sorted_clusters.csv; fi",
            rf"awk 'NR>1 {{print}}' {assembly_dir}/{OUTPUT_DIRS[7]}/{basename}/{basename}_summary.csv | sort -t ',' -k1,1 > {assembly_dir}/{OUTPUT_DIRS[9]}/{basename}/{basename}_sorted_summary.csv",
            rf"""
            join -t',' -1 1 -2 1 {assembly_dir}/{OUTPUT_DIRS[9]}/{basename}/{basename}_sorted_clusters.csv \
            {assembly_dir}/{OUTPUT_DIRS[9]}/{basename}/{basename}_sorted_summary.csv | \
            awk -F',' 'BEGIN {{print "contig_id,cluster_id,contig_length,gene_count,viral_genes,host_genes,checkv_quality,miuvig_quality,completeness,completeness_method,contamination,kmer_freq,score"}} \
            {{contig_length=$3; gene_count=$4; viral_genes=$5; host_genes=$6; miuvig_quality=$8; completeness=$9; contamination=$11; kmer_freq=$12; score=kmer_freq*viral_genes; print $0","score}}' > \
            {assembly_dir}/{OUTPUT_DIRS[9]}/{basename}/{basename}_final_summary.csv
            """,
            rf"awk -F',' 'NR==1 {{print; next}} {{if($13 > max[$2]) {{max[$2]=$13; best[$2]=$0}}}} END {{for(c in best) print best[c]}}' {assembly_dir}/{OUTPUT_DIRS[9]}/{basename}/{basename}_final_summary.csv > {assembly_dir}/{OUTPUT_DIRS[9]}/{basename}/{basename}_selected_genomes.csv",
            rf"if [ ! -f {assembly_dir}/{OUTPUT_DIRS[9]}/{basename}/{basename}_selected_genomes.fasta ]; then awk -F',' 'NR>1 {{print $1}}' {assembly_dir}/{OUTPUT_DIRS[9]}/{basename}/{basename}_selected_genomes.csv | while read contig_id; do cat {assembly_dir}/{OUTPUT_DIRS[8]}/{basename}/$contig_id.fa >> {assembly_dir}/{OUTPUT_DIRS[9]}/{basename}/{basename}_selected_genomes.fasta; done; fi",
            rf"rm {assembly_dir}/{OUTPUT_DIRS[9]}/{basename}/{basename}_sorted_clusters.csv {assembly_dir}/{OUTPUT_DIRS[9]}/{basename}/{basename}_sorted_summary.csv"
        ],
        10: [ # Depth and coverage
            f"echo \"Tool version (samtools):\" $(samtools version | head -n 1)",
            f"echo \"Tool version (minimap2):\" $(minimap2 --version)",
            f"minimap2 -ax map-ont {assembly_dir}/{OUTPUT_DIRS[9]}/{basename}/{basename}_selected_genomes.fasta {assembly_dir}/filtered_{basename}.fastq | samtools sort -o {assembly_dir}/{OUTPUT_DIRS[10]}/{basename}.bam",
            f"samtools index {assembly_dir}/{OUTPUT_DIRS[10]}/{basename}.bam -o {assembly_dir}/{OUTPUT_DIRS[10]}/{basename}.bam.bai",
            f"samtools coverage {assembly_dir}/{OUTPUT_DIRS[10]}/{basename}.bam > {assembly_dir}/{OUTPUT_DIRS[10]}/{basename}_aln_temp.tsv",
            f"awk 'BEGIN{{FS=OFS=\"\\t\"}} NR==1{{$1=\"contig\"; for(i=2;i<=NF;i++) $i=\"samtools_\"$i}} 1' {assembly_dir}/{OUTPUT_DIRS[10]}/{basename}_aln_temp.tsv > {assembly_dir}/{OUTPUT_DIRS[10]}/{basename}_aln.tsv"
        ],
        11: [ # Cleanup temporary files
            rf"rm {assembly_dir}/filtered_{basename}.fastq {assembly_dir}/{OUTPUT_DIRS[10]}/{basename}_aln_temp.tsv {assembly_dir}/{OUTPUT_DIRS[10]}/{basename}.bam {assembly_dir}/{OUTPUT_DIRS[10]}/{basename}.bam.bai"
        ]
    }
    return commands_by_step.get(step_number, [])

def run_assembly(output_dir, input_files, task_params):
    assembly_dir = os.path.join(output_dir, "assembly_task")

    # Ensure the assembly directory and subdirectories exist
    os.makedirs(assembly_dir, exist_ok=True)
    log_path = os.path.join(output_dir, "assembly.log")
    setup_logging(log_path)

    logger.info(f"Using folder: \"{assembly_dir}\"")

    for sub_dir in OUTPUT_DIRS.values():
        os.makedirs(os.path.join(assembly_dir, sub_dir), exist_ok=True)

    # Start total timer
    total_start_time = time.time()

    # Define used environments per step
    envs = {
        0: {"env_name": "QPPL_env"},
        1: {"env_name": "QPPL_env"},
        2: {"env_name": "QPPL_env2"},
        3: {"env_name": "QPPL_env"},
        4: {"env_name": "QPPL_env"},
        5: {"env_name": "QPPL_env"},
        6: {"env_name": "QPPL_env2"},
        7: {"env_name": "QPPL_env"},
        8: {"env_name": "QPPL_env"},
        9: {"env_name": "QPPL_env4"},
        10: {"env_name": "QPPL_env"},
        11: {"env_name": "QPPL_env"}
    }

    # Process each input file
    for file in input_files:
        basename = re.sub(r'(\.fastq|\.fastq\.gz|\.fq|\.fq\.gz)$', '', file, flags=re.IGNORECASE)
        
        # Start task timer per file
        start_time = time.time()

        logger.info("=" * 60)
        logger.info(f"Starting assembly tasks for file: {file}")

        # Run commands in environments
        for step_number in sorted(envs.keys()):
            commands = generate_commands(assembly_dir, output_dir, file, task_params, basename, step_number)
            for command in commands:
                run_command(envs[step_number]['env_name'], command, file, step_number)

        # Log time taken for the file
        elapsed = time.time() - start_time
        logger.info(f"Finished assembly tasks for file: {file} (elapsed time: {elapsed:.2f} seconds)")
    
    # Calculate total elapsed time
    total_elapsed = time.time() - total_start_time
    logger.info(f"Total elapsed time for all tasks: {total_elapsed:.2f} seconds")