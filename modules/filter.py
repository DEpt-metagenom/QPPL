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

            
def run_filter(input_dir, output_dir, prefix, input_files, task_params):

    filter_dir = os.path.join(output_dir, "filter_task")

    if os.path.isdir(filter_dir):
        print(f"\033[1;31m[INFO]\033[0m Detected existing folder: \"{filter_dir}\"")
    else:
        print(f"\033[1;32m[INFO]\033[0m Created new folder: \"{filter_dir}\"")
        os.makedirs(filter_dir, exist_ok=True)

    ### LOOP OVER ALL INPUT FILES
    for file in input_files:
        basename = re.sub(r'(\.fastq|\.fastq\.gz|\.fq|\.fq\.gz)$', '', file)

        output_file_porechop = os.path.join(filter_dir, f"noadapter_{basename}.fastq")
        output_file_nanotools = os.path.join(filter_dir, f"filtered_{file}")

        ### DEBUG
        print(f"\033[1;34m[DEBUG]\033[0m input_dir: {input_dir}")
        print(f"\033[1;34m[DEBUG]\033[0m output_dir: {output_dir}")
        print(f"\033[1;34m[DEBUG]\033[0m prefix: {prefix}")
        print(f"\033[1;34m[DEBUG]\033[0m file: {file}")
        print(f"\033[1;34m[DEBUG]\033[0m basename: {basename}")
        print(f"\033[1;34m[DEBUG]\033[0m porechop_out: {output_file_porechop}")
        print(f"\033[1;34m[DEBUG]\033[0m nanotools_out: {output_file_nanotools}")

        for key, value in task_params.get("filter", {}).items():
            print(f"\033[1;34m[DEBUG]\033[0m {key}: {value}")

        ### RUN FILTERING TOOLS
        if task_params['filter']['cut_adapter'] == 'y':
            envs = {
                1: {
                    "env_name": "QPPL_env",
                    "commands": [
                        f"porechop_abi -abi --discard_middle -i {input_dir}/{file} -o {output_file_porechop}",
                        f"gzip {output_file_porechop}"
                    ]
                },
                2: {
                    "env_name": "QPPL_env",
                    "commands": [
                        f"if [ ! -f {filter_dir}/lambda.fasta ]; then wget -q https://github.com/wdecoster/nanolyse/raw/master/reference/lambda.fasta.gz -O {filter_dir}/lambda.fasta.gz && gunzip {filter_dir}/lambda.fasta.gz; fi",
                        f"if [ ! -f {filter_dir}/DNA_CS.fasta ]; then wget -q https://raw.githubusercontent.com/wdecoster/nanolyse/master/reference/DNA_CS.fasta -O {filter_dir}/DNA_CS.fasta; fi",
                        f"gunzip -c {output_file_porechop} | NanoLyse -r {filter_dir}/lambda.fasta | NanoLyse -r {filter_dir}/DNA_CS.fasta | NanoFilt -q {task_params['filter']['filter_quality']} -l {task_params['filter']['filter_length']} | gzip > {output_file_nanotools}",
                        f"rm NanoLyse.log",
                        f"rm -r tmp"
                    ]
                }
            }
        else:
            envs = {
                1: {
                    "env_name": "QPPL_env",
                    "commands": [
                        f"if [ ! -f {filter_dir}/lambda.fasta ]; then wget -q https://github.com/wdecoster/nanolyse/raw/master/reference/lambda.fasta.gz -O {filter_dir}/lambda.fasta.gz && gunzip {filter_dir}/lambda.fasta.gz; fi",
                        f"if [ ! -f {filter_dir}/DNA_CS.fasta ]; then wget -q https://raw.githubusercontent.com/wdecoster/nanolyse/master/reference/DNA_CS.fasta -O {filter_dir}/DNA_CS.fasta; fi",
                        f"gunzip -c {input_dir}/{file} | NanoLyse -r {filter_dir}/lambda.fasta | NanoLyse -r {filter_dir}/DNA_CS.fasta | NanoFilt -q {task_params['filter']['filter_quality']} -l {task_params['filter']['filter_length']} | gzip > {output_file_nanotools}",
                        f"rm NanoLyse.log"
                    ]
                }
            }


        for step_number in sorted(envs.keys()):
            env = envs[step_number]
            run_commands_in_env(env["env_name"], env["commands"], file, step_number)