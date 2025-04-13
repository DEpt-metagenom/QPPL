import os
import subprocess
import re
import modules.parser
import modules.version
import modules.config
import modules.filter
import modules.assembly
import modules.characterization
import modules.readquality

def main():
    parser = modules.parser.create_parser()
    args = parser.parse_args()

    if args.version:
        print(modules.version.LOGO.strip())
    elif args.generate_config:
        modules.config.generate_default_config()
    elif args.config:
        try:
            input_dir, output_dir, prefix, task, task_params = modules.config.read_config(args.config)

            # Print the configuration
            print(f"Input Directory: {input_dir}")
            print(f"Output Directory: {output_dir}")
            print(f"Prefix: {prefix}")
            print(f"Task: {task}")
            print(f"Task Parameters: {task_params}")

            if task == 'readquality':
                input_files = [
                    f for f in os.listdir(input_dir) 
                    if f.endswith((".fastq", ".fastq.gz", ".fq", ".fq.gz"))
                ]

                if not input_files:
                    print("\033[1;31m[ERROR]\033[0m No FASTQ files found in input directory.")
                    return

                modules.readquality.run_read_quality(input_dir, output_dir, prefix, input_files, task_params)

            elif task == 'filter':
                input_files = [
                    f for f in os.listdir(input_dir) 
                    if f.endswith((".fastq", ".fastq.gz", ".fq", ".fq.gz"))
                ]

                if not input_files:
                    print("\033[1;31m[ERROR]\033[0m No FASTQ files found in input directory.")
                    return

                modules.filter.run_filter(input_dir, output_dir, prefix, input_files, task_params)

            elif task == 'assembly':
                input_files = [
                    f for f in os.listdir(input_dir) 
                    if f.endswith((".fastq", ".fastq.gz", ".fq", ".fq.gz"))
                ]

                if not input_files:
                    print("\033[1;31m[ERROR]\033[0m No FASTQ files found in input directory.")
                    return

                modules.assembly.run_assembly(input_dir, output_dir, prefix, input_files, task_params)
            
            elif task == 'characterization':
                input_files = [
                    f for f in os.listdir(input_dir) 
                    if f.endswith((".fastq", ".fastq.gz", ".fq", ".fq.gz"))
                ]

                if not input_files:
                    print("\033[1;31m[ERROR]\033[0m No FASTQ files found in input directory.")
                    return

                modules.characterization.run_characterization(input_dir, output_dir, prefix, input_files, task_params)
            
            elif task == 'all':
                input_files = [
                    f for f in os.listdir(input_dir) 
                    if f.endswith((".fastq", ".fastq.gz", ".fq", ".fq.gz"))
                ]

                if not input_files:
                    print("\033[1;31m[ERROR]\033[0m No FASTQ files found in input directory.")
                    return

                modules.filter.run_filter(input_dir, output_dir, prefix, input_files, task_params)
                modules.assembly.run_assembly(input_dir, output_dir, prefix, input_files, task_params)
                modules.characterization.run_characterization(input_dir, output_dir, prefix, input_files, task_params)

        except ValueError as e:
            print(f"Error reading config: {e}")

if __name__ == "__main__":
    main()