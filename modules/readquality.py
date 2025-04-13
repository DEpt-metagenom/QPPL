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

def run_read_quality(input_dir, output_dir, prefix, input_files, task_params):

    read_quality_dir = os.path.join(output_dir, "read_quality_task")

    ### CREATE READ QUALITY DIRECTORY
    if os.path.isdir(read_quality_dir):
        print(f"\033[1;31m[INFO]\033[0m Detected existing folder: \"{read_quality_dir}\"")
    else:
        print(f"\033[1;32m[INFO]\033[0m Created new folder: \"{read_quality_dir}\"")
        os.makedirs(read_quality_dir, exist_ok=True)

    ### LOOP OVER ALL INPUT FILES
    for file in input_files:
        basename = re.sub(r'(\.fastq|\.fastq\.gz|\.fq|\.fq\.gz)$', '', file)

        ### PRINT DEBUG INFORMATION
        print(f"\033[1;34m[DEBUG]\033[0m input_dir: {input_dir}")
        print(f"\033[1;34m[DEBUG]\033[0m output_dir: {output_dir}")
        print(f"\033[1;34m[DEBUG]\033[0m prefix: {prefix}")
        print(f"\033[1;34m[DEBUG]\033[0m file: {file}")
        print(f"\033[1;34m[DEBUG]\033[0m basename: {basename}")

        ### RUN READ QUALITY TOOLS
        envs = {
            1: {
                "env_name": "QPPL_env",
                "commands": [
                    f"nanoQC -o {read_quality_dir}/ {input_dir}/{file}",
                    f"mv {read_quality_dir}/NanoQC.log {read_quality_dir}/{basename}_NanoQC.log",
                    f"mv {read_quality_dir}/nanoQC.html {read_quality_dir}/{basename}_nanoQC.html",
                    f"nanoq -i {input_dir}/{file} -s -H > {read_quality_dir}/{basename}_nanoq.log"
                ]
            }
        }

        for step_number in sorted(envs.keys()):
            env = envs[step_number]
            run_commands_in_env(env["env_name"], env["commands"], file, step_number)