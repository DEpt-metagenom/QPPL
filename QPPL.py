import os
import sys
import logging
import modules.parser
import modules.version
import modules.config
import modules.filter
import modules.assembly
import modules.characterization
import modules.readquality
import modules.setup
import modules.databases

# Logger setup with color support
class ColorFormatter(logging.Formatter):
    COLORS = {
        logging.DEBUG: "\033[36m",    # Cyan
        logging.INFO: "\033[32m",     # Green
        logging.WARNING: "\033[33m",  # Yellow
        logging.ERROR: "\033[31m",    # Red
        logging.CRITICAL: "\033[41m", # Red background
    }
    RESET = "\033[0m"

    def format(self, record):
        color = self.COLORS.get(record.levelno, self.RESET)
        message = super().format(record)
        return f"{color}{message}{self.RESET}" if sys.stdout.isatty() else message

# Logger configuration
logger = logging.getLogger("QPPL")
logger.setLevel(logging.DEBUG)
console_handler = logging.StreamHandler()
console_handler.setFormatter(ColorFormatter("[%(levelname)s] %(asctime)s - %(message)s", "%Y-%m-%d %H:%M:%S"))
logger.addHandler(console_handler)

# Function to get FASTQ files from input directory
def get_fastq_files(input_dir):
    return [f for f in os.listdir(input_dir) if f.endswith((".fastq", ".fastq.gz", ".fq", ".fq.gz"))]

# Main function
def main():
    parser = modules.parser.create_parser()
    args = parser.parse_args()

    if args.version:
        print(modules.version.LOGO.strip())
        return
    
    if args.generate_config:
        modules.config.generate_default_config()
        return

    if args.install_envs:
        modules.setup.run_install_envs(os.path.expanduser(args.vhulk_dir))
        return

    if args.download_databases:
        modules.databases.run_download_databases(args.databases_dir)
        return

    try:
        input_dir, output_dir, task, task_params = modules.config.read_config(args.config)

        # Print the configuration
        print(f"Input Directory: {input_dir}\nOutput Directory: {output_dir}\nTask: {task}\nTask Parameters: {task_params}")

        input_files = get_fastq_files(input_dir)
        if not input_files:
            logger.error("No FASTQ files found in input directory.")
            return

        # Execute the task or all tasks
        if task == 'readquality':
            modules.readquality.run_read_quality(input_dir, output_dir, input_files, task_params)
        elif task == 'filter':
            modules.filter.run_filter(input_dir, output_dir, input_files, task_params)
        elif task == 'assembly':
            modules.assembly.run_assembly(output_dir, input_files, task_params)
        elif task == 'characterization':
            modules.characterization.run_characterization(output_dir, input_files, task_params)
        elif task == 'all':
            # Execute tasks in the specified order
            modules.filter.run_filter(input_dir, output_dir, input_files, task_params)
            modules.assembly.run_assembly(output_dir, input_files, task_params)
            modules.characterization.run_characterization(output_dir, input_files, task_params)
        else:
            logger.error(f"Unknown task: {task}")

    except ValueError as e:
        logger.error(f"Error reading config: {e}")

if __name__ == "__main__":
    main()