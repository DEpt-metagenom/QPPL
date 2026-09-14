#!/usr/bin/env python3
"""
run.py - runs QPPL.py under SLURM via srun.

Just run:

    python3 run.py

and it requests the allocation and launches QPPL.py itself
(equivalent to running the srun command below by hand).
"""
import os
import subprocess
import sys

QPPL_DIR = os.path.dirname(os.path.abspath(__file__))

SRUN_CMD = [
    "srun",
    "--nodes=1",
    "--cpus-per-task=4",  # TEMP: lowered from 10 to fit local laptop test (4 CPUs avail)
    "--mem=24G",  # TEMP: lowered from 120G to fit local laptop test (~30G avail)
    "python3", "QPPL.py", "--config", "qppl.conf",
]

try:
    subprocess.run(SRUN_CMD, check=True, cwd=QPPL_DIR)
except FileNotFoundError:
    sys.exit("ERROR: 'srun' not found on PATH -- run this from a machine with SLURM client tools available.")
