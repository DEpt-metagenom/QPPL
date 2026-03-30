import argparse
from .version import LOGO

class CustomArgumentParser(argparse.ArgumentParser):
    def format_help(self):
        default_help = super().format_help()
        return f"{LOGO.strip()}\n\n{default_help}"

def create_parser():
    parser = CustomArgumentParser(description="QPPL (Quick Phage PipeLine) help menu:", add_help=False)

    parser.add_argument(
        "-h", "--help",
        action="help",
        help="Show this help message and exit"
    )
    parser.add_argument(
        "-v", "--version",
        action="store_true",
        help="Display QPPL version and logo"
    )
    parser.add_argument(
        "-gc", "--generate-config",
        action="store_true",
        help="Generate a default configuration file (overrides existing)"
    )
    parser.add_argument(
        "-c", "--config",
        type=str,
        default='qppl.conf',
        help="Path to the configuration file"
    )
    return parser