import os
import subprocess
import re

def run_commands_in_env(env_name, commands, file, env_number):
    for command in commands:
        print(f"\n\033[1;32m[INFO] Running command in {env_name} on {file} (step {env_number}):\n  \033[0m$ {command}")
        result = subprocess.run(f"conda run -n {env_name} bash -c \"{command}\"", 
                                shell=True, text=True, capture_output=True)
        if result.returncode != 0:
            print(f"\033[1;31m[ERROR] Command failed: {command}\n\033[0m{result.stderr}")
        else:
            print(result.stdout)

def run_assembly(input_dir, output_dir, prefix, input_files, task_params):

    assembly_dir = os.path.join(output_dir, "assembly_task")

    sub_dirs = [
        "out_flye", "out_goldrush", "out_raven", 
        "out_shasta", "out_wtdbg2", "out_polish", 
        "out_checkv", "contigs", "out_vclust"
    ]

    ### CREATE ASSEMBLY DIRECTORY AND SUBDIRECTORIES
    if os.path.isdir(assembly_dir):
        print(f"\033[1;31m[INFO]\033[0m Detected existing folder: \"{assembly_dir}\"")
    else:
        print(f"\033[1;32m[INFO]\033[0m Created new folder: \"{assembly_dir}\"")
        os.makedirs(assembly_dir, exist_ok=True)

    for sub_dir in sub_dirs:
        path = os.path.join(assembly_dir, sub_dir)
        if not os.path.isdir(path):
            print(f"\033[1;32m[INFO]\033[0m Created new folder: \"{path}\"")
            os.makedirs(path, exist_ok=True)
        else:
            print(f"\033[1;31m[INFO]\033[0m Folder already exists: \"{path}\"")

    ### LOOP OVER ALL INPUT FILES
    for file in input_files:
        basename = re.sub(r'(\.fastq|\.fastq\.gz|\.fq|\.fq\.gz)$', '', file)

        ### PRINT DEBUG INFORMATION
        print(f"\033[1;34m[DEBUG]\033[0m input_dir: {input_dir}")
        print(f"\033[1;34m[DEBUG]\033[0m output_dir: {output_dir}")
        print(f"\033[1;34m[DEBUG]\033[0m prefix: {prefix}")
        print(f"\033[1;34m[DEBUG]\033[0m file: {file}")
        print(f"\033[1;34m[DEBUG]\033[0m basename: {basename}")

        for key, value in task_params.get("assembly", {}).items():
            print(f"\033[1;34m[DEBUG]\033[0m {key}: {value}")

        ### RUN ASSEMBLY TOOLS
        meta_flag = "--meta" if task_params['assembly']['flye_meta'] == 'y' else ""
        wtdbg2_realign_flag = "-R" if task_params['assembly']['wtdbg2_realign'] == 'y' else ""

        envs = {
            1: {
                "env_name": "QPPL_env",
                "commands": [
                    f"flye --{task_params['assembly']['flye_mode']} {output_dir}/filter_task/filtered_{file} -o {assembly_dir}/out_flye/{basename} --threads {task_params['assembly']['flye_threads']} {meta_flag}",
                    f"seqkit replace -p '.+' -r flye_contig_{{nr}} {assembly_dir}/out_flye/{basename}/assembly.fasta | seqkit seq -w 0  > {assembly_dir}/out_flye/{basename}_flye.fasta",
                    f"rm -r {assembly_dir}/out_flye/{basename}",

                    f"raven -t {task_params['assembly']['raven_threads']} {output_dir}/filter_task/filtered_{file} > {assembly_dir}/out_raven/{basename}_temp.fasta",
                    f"seqkit replace -p '.+' -r raven_contig_{{nr}} {assembly_dir}/out_raven/{basename}_temp.fasta | seqkit seq -w 0 > {assembly_dir}/out_raven/{basename}_raven.fasta",
                    f"rm {assembly_dir}/out_raven/{basename}_temp.fasta",
                    f"rm raven.cereal",

                    f"mkdir {assembly_dir}/out_wtdbg2/{basename}",
                    f"wtdbg2 -g {task_params['assembly']['wtdbg2_genome_size']} -t {task_params['assembly']['wtdbg2_threads']} -i {output_dir}/filter_task/filtered_{file} -fo {assembly_dir}/out_wtdbg2/{basename}/{basename} -L {task_params['assembly']['wtdbg2_read_length']} -k {task_params['assembly']['wtdbg2_kmer_fsize']} -p {task_params['assembly']['wtdbg2_kmer_psize']} -s {task_params['assembly']['wtdbg2_min_similarity']} {wtdbg2_realign_flag}",
                    f"wtpoa-cns -t {task_params['assembly']['wtdbg2_threads']} -i {assembly_dir}/out_wtdbg2/{basename}/{basename}.ctg.lay.gz -fo {assembly_dir}/out_wtdbg2/{basename}/{basename}.ctg.fasta",
                    f"seqkit replace -p '.+' -r wtdbg2_contig_{{nr}} {assembly_dir}/out_wtdbg2/{basename}/{basename}.ctg.fasta | seqkit seq -w 0 > {assembly_dir}/out_wtdbg2/{basename}_wtdbg2.fasta",
                    f"rm -r {assembly_dir}/out_wtdbg2/{basename}",

                    f"if [ ! -f {assembly_dir}/Nanopore-May2022.conf ]; then wget -q https://raw.githubusercontent.com/paoloshasta/shasta/main/conf/Nanopore-May2022.conf -O {assembly_dir}/Nanopore-May2022.conf; fi",
                    f"sed -i 's/guppy-5.0.7-b/guppy-5.0.7-a/g' {assembly_dir}/Nanopore-May2022.conf",
                    f"sed -i 's/minReadLength = 10000/minReadLength = 1000/g' {assembly_dir}/Nanopore-May2022.conf",
                    f"gunzip -c {output_dir}/filter_task/filtered_{file} > {assembly_dir}/out_shasta/{basename}.fq",
                    f"seqkit grep -v -s -r -p '[^ACGT]' {assembly_dir}/out_shasta/{basename}.fq > {assembly_dir}/out_shasta/{basename}_clean.fq",
                    f"shasta --input {assembly_dir}/out_shasta/{basename}_clean.fq --config {assembly_dir}/Nanopore-May2022.conf --threads {task_params['assembly']['shasta_threads']} --assemblyDirectory {assembly_dir}/out_shasta/{basename}",
                    f"seqkit replace -p '.+' -r shasta_contig_{{nr}} {assembly_dir}/out_shasta/{basename}/Assembly.fasta | seqkit seq -w 0 > {assembly_dir}/out_shasta/{basename}_shasta.fasta",
                    f"rm {assembly_dir}/out_shasta/{basename}.fq",
                    f"rm {assembly_dir}/out_shasta/{basename}_clean.fq",
                    f"rm -r {assembly_dir}/out_shasta/{basename}"
                ]
            },
            2: {
                "env_name": "QPPL_env2",
                "commands": [
                    f"mkdir {assembly_dir}/out_goldrush/{basename}",
                    f"gunzip -c {output_dir}/filter_task/filtered_{file} > {assembly_dir}/out_goldrush/{basename}/tmp.fq",
                    f"cd {assembly_dir}/out_goldrush/{basename} && goldrush run reads=tmp p=tmp G={task_params['assembly']['goldrush_genome_size']} m={task_params['assembly']['goldrush_min_read_length']} t={task_params['assembly']['goldrush_threads']} && cd -",
                    f"mv {assembly_dir}/out_goldrush/{basename}/goldrush_intermediate_files/tmp*.polished.fa {assembly_dir}/out_goldrush/{basename}/assembly.fasta",
                    f"seqkit replace -p '.+' -r goldrush_contig_{{nr}} {assembly_dir}/out_goldrush/{basename}/assembly.fasta | seqkit seq -w 0 > {assembly_dir}/out_goldrush/{basename}_goldrush.fasta",
                    f"rm -r {assembly_dir}/out_goldrush/{basename}",
                ]
            },
            3: {
                "env_name": "QPPL_env2",
                "commands": [
                    rf"gunzip -c {output_dir}/filter_task/filtered_{file} > {assembly_dir}/out_polish/tmp_{basename}.fastq",
                    rf"for asm in flye raven wtdbg2 shasta goldrush; do minimap2 -x map-ont -t 8 -a {assembly_dir}/out_\$asm/{basename}_\$asm.fasta {assembly_dir}/out_polish/tmp_{basename}.fastq > {assembly_dir}/out_polish/{basename}_align_\$asm.sam; done",
                    rf"for asm in flye raven wtdbg2 shasta goldrush; do racon -t 8 -f {assembly_dir}/out_polish/tmp_{basename}.fastq {assembly_dir}/out_polish/{basename}_align_\$asm.sam {assembly_dir}/out_\$asm/{basename}_\$asm.fasta > {assembly_dir}/out_polish/{basename}_asmbl_\$asm.fa; done",
                    rf"for asm in flye raven wtdbg2 shasta goldrush; do mini_align -i {assembly_dir}/out_polish/tmp_{basename}.fastq -r {assembly_dir}/out_polish/{basename}_asmbl_\$asm.fa -P -m -p {assembly_dir}/out_polish/{basename}_align_\$asm -t 8; done",
                    rf"sleep 5",
                    rf"for asm in flye raven wtdbg2 shasta goldrush; do medaka inference --model r1041_e82_400bps_sup_g615 {assembly_dir}/out_polish/{basename}_align_\$asm.bam {assembly_dir}/out_polish/{basename}_align_\$asm.hdf; done",
                    rf"sleep 5",
                    rf"for asm in flye raven wtdbg2 shasta goldrush; do medaka sequence {assembly_dir}/out_polish/{basename}_align_\$asm.hdf {assembly_dir}/out_polish/{basename}_asmbl_\$asm.fa {assembly_dir}/out_polish/{basename}_polished_\$asm.fa; done",
                    rf"find {assembly_dir}/out_polish/ -type f ! -name \"*_polished_*.fa\" -exec rm {{}} \;",
                ]
            },
            4: {
                "env_name": "CheckV",
                "commands": [
                    rf"mkdir {assembly_dir}/out_checkv/{basename}",
                    rf"for asm in flye shasta raven wtdbg2 goldrush; do cat {assembly_dir}/out_polish/{basename}_polished_\$asm.fa >> {assembly_dir}/out_checkv/{basename}/{basename}_polished.fa; done",
                    rf"checkv end_to_end {assembly_dir}/out_checkv/{basename}/{basename}_polished.fa {assembly_dir}/out_checkv/{basename} -d {task_params['assembly']['checkv_db']} -t {task_params['assembly']['checkv_threads']}"
                ]
            },
            5: {
                "env_name": "QPPL_env",
                "commands": [
                    rf"csvcut -t -c 1,2,5,6,7,8,9,10,11,12,13 {assembly_dir}/out_checkv/{basename}/quality_summary.tsv | head -n 1 > {assembly_dir}/out_checkv/{basename}/{basename}_summary.csv",
                    rf"csvcut -t -c 1,2,5,6,7,8,9,10,11,12,13 {assembly_dir}/out_checkv/{basename}/quality_summary.tsv | csvgrep -c 7 -m High-quality | tail -n +2 >> {assembly_dir}/out_checkv/{basename}/{basename}_summary.csv",
                    rf"mkdir {assembly_dir}/contigs/{basename}",
                    rf"csvcut -c 1 {assembly_dir}/out_checkv/{basename}/{basename}_summary.csv | tail -n +2 | while IFS= read -r contig; do seqtk subseq {assembly_dir}/out_checkv/{basename}/{basename}_polished.fa <(echo \"\$contig\") > {assembly_dir}/contigs/{basename}/\$contig.fa; done",
                ]
            },
            6: {
                "env_name": "QPPL_env4",
                "commands": [
                    rf"mkdir {assembly_dir}/out_vclust/{basename}",
                    rf"cat {assembly_dir}/contigs/{basename}/*.fa > {assembly_dir}/out_vclust/{basename}/{basename}_contigs.fa",
                    rf"if [ \$(grep -c \">\" {assembly_dir}/out_vclust/{basename}/{basename}_contigs.fa) -eq 1 ]; then cp {assembly_dir}/out_vclust/{basename}/{basename}_contigs.fa {assembly_dir}/out_vclust/{basename}/{basename}_selected_genomes.fasta; fi",
                    rf"if [ ! -f {assembly_dir}/out_vclust/{basename}/{basename}_selected_genomes.fasta ]; then vclust prefilter -i {assembly_dir}/out_vclust/{basename}/{basename}_contigs.fa -o {assembly_dir}/out_vclust/{basename}/{basename}_fltr.txt --min-ident 0.7; fi",
                    rf"if [ ! -f {assembly_dir}/out_vclust/{basename}/{basename}_selected_genomes.fasta ]; then vclust align -i {assembly_dir}/out_vclust/{basename}/{basename}_contigs.fa -o {assembly_dir}/out_vclust/{basename}/{basename}_align.tsv --filter {assembly_dir}/out_vclust/{basename}/{basename}_fltr.txt --outfmt lite --out-ani 0.95 --out-qcov 0.85; fi",
                    rf"if [ ! -f {assembly_dir}/out_vclust/{basename}/{basename}_selected_genomes.fasta ]; then vclust cluster -i {assembly_dir}/out_vclust/{basename}/{basename}_align.tsv -o {assembly_dir}/out_vclust/{basename}/{basename}_clusters.tsv --ids {assembly_dir}/out_vclust/{basename}/{basename}_align.ids.tsv --algorithm complete --metric ani --ani 0.95 --qcov 0.85; fi",
                    rf"if [ ! -f {assembly_dir}/out_vclust/{basename}/{basename}_selected_genomes.fasta ]; then tail -n +2 {assembly_dir}/out_vclust/{basename}/{basename}_clusters.tsv | sort -k1,1 | tr \"\t\" \",\" > {assembly_dir}/out_vclust/{basename}/{basename}_sorted_clusters.csv; else contig=\$(head -n 1 {assembly_dir}/out_vclust/{basename}/{basename}_selected_genomes.fasta | sed 's/>//g'); echo \"\${{contig}},0\" > {assembly_dir}/out_vclust/{basename}/{basename}_sorted_clusters.csv; fi",
                    rf"tail -n +2 {assembly_dir}/out_checkv/{basename}/{basename}_summary.csv | sort -t ',' -k1,1 > {assembly_dir}/out_vclust/{basename}/{basename}_sorted_summary.csv",
                    rf"join -t',' -1 1 -2 1 {assembly_dir}/out_vclust/{basename}/{basename}_sorted_clusters.csv {assembly_dir}/out_vclust/{basename}/{basename}_sorted_summary.csv > {assembly_dir}/out_vclust/{basename}/{basename}_tmp.csv",
                    rf"echo \"contig_id,cluster_id,contig_length,gene_count,viral_genes,host_genes,checkv_quality,miuvig_quality,completeness,completeness_method,contamination,kmer_freq\" > {assembly_dir}/out_vclust/{basename}/{basename}_tmp2.csv",
                    rf"cat {assembly_dir}/out_vclust/{basename}/{basename}_tmp.csv >> {assembly_dir}/out_vclust/{basename}/{basename}_tmp2.csv",
                    rf"awk -F',' 'NR==1 {{print \$0\",score\"; next}} {{print \$0\",\"((\$3/\$12)/\$5)*(\$9/100)*((100-\$11)/100)}}' {assembly_dir}/out_vclust/{basename}/{basename}_tmp2.csv > {assembly_dir}/out_vclust/{basename}/{basename}_final_summary.csv",
                    rf"awk -F',' 'NR==1 {{print \$0; next}} {{if(\$13 > max[\$2]) {{max[\$2]=\$13; best[\$2]=\$0}}}} END {{for(c in best) print best[c]}}' {assembly_dir}/out_vclust/{basename}/{basename}_final_summary.csv > {assembly_dir}/out_vclust/{basename}/{basename}_selected_genomes.csv",
                    rf"if [ ! -f {assembly_dir}/out_vclust/{basename}/{basename}_selected_genomes.fasta ]; then awk -F',' 'NR > 1 {{print \$1}}' {assembly_dir}/out_vclust/{basename}/{basename}_selected_genomes.csv | while read contig_id; do cat {assembly_dir}/contigs/{basename}/\$contig_id.fa >> {assembly_dir}/out_vclust/{basename}/{basename}_selected_genomes.fasta; done; fi",
                    rf"rm {assembly_dir}/out_vclust/{basename}/{basename}_tmp.csv",
                    rf"rm {assembly_dir}/out_vclust/{basename}/{basename}_tmp2.csv",
                ]
            }
        }

        for step_number in sorted(envs.keys()):
            env = envs[step_number]
            run_commands_in_env(env["env_name"], env["commands"], file, step_number)
