import os
import re
import struct
import subprocess
import logging
import time

# Call logger
logger = logging.getLogger("QPPL")

VHULK_REPO_URL = "https://github.com/LaboratorioBioinformatica/vHULK.git"

# (env_name, yml file relative to repo root)
ENV_SPECS = [
    ("QPPL_env", "envs/QPPL_env.yml"),
    ("QPPL_env2", "envs/QPPL_env2.yml"),
    ("QPPL_env3", "envs/QPPL_env3.yml"),
    ("QPPL_env4", "envs/QPPL_env4.yml"),
    ("vHULK", "envs/vHULK_env.yml"),
]

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

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

# Function to run a command with live (inherited) stdout/stderr, since env
# creation and package installs can take many minutes and shouldn't look hung.
def run_live(command, description):
    log_to_file_only(f"Running: $ {' '.join(command)}")
    logger.info(description)
    result = subprocess.run(command)
    if result.returncode != 0:
        logger.error(f"Failed ({result.returncode}): {' '.join(command)}")
        return False
    return True

# Function to capture a command's output quietly, for version/state checks
def run_quiet(command):
    result = subprocess.run(command, text=True, capture_output=True)
    return result.returncode, result.stdout.strip(), result.stderr.strip()

# Function to check whether a conda environment already exists
def conda_env_exists(env_name):
    returncode, stdout, _ = run_quiet(["conda", "env", "list"])
    if returncode != 0:
        return False
    names = [line.split()[0] for line in stdout.splitlines() if line and not line.startswith("#")]
    return env_name in names

# Function to create a conda environment from a yml file if it doesn't already exist
def create_env(env_name, yml_relpath):
    if conda_env_exists(env_name):
        logger.info(f"Conda environment '{env_name}' already exists, skipping creation.")
        return True

    yml_path = os.path.join(REPO_ROOT, yml_relpath)
    logger.info(f"Creating conda environment '{env_name}' from {yml_path} (this can take a while)...")
    return run_live(["conda", "env", "create", "-n", env_name, "-f", yml_path], f"conda env create -n {env_name}")

# Function to get an installed pip package's version inside a conda env, or None
def _pip_version(env_name, package):
    returncode, stdout, _ = run_quiet(["conda", "run", "-n", env_name, "pip", "show", package])
    if returncode != 0:
        return None
    match = re.search(r"^Version:\s*(\S+)", stdout, re.MULTILINE)
    return match.group(1) if match else None

# Function to parse a simple dotted version string into a comparable tuple
def _version_tuple(version):
    parts = re.findall(r"\d+", version or "")
    return tuple(int(p) for p in parts) if parts else (0,)

# pharokka 1.7.3 pins mmseqs2<14 and dnaapler 0.8.1, incompatible with pharokka
# 1.8.0's requirements. Bump all three in one guarded pass.
def fix_qppl_env():
    env_name = "QPPL_env"
    if not conda_env_exists(env_name):
        logger.warning(f"Conda environment '{env_name}' does not exist, skipping its fixes.")
        return

    pharokka_version = _pip_version(env_name, "pharokka")
    if pharokka_version and _version_tuple(pharokka_version) >= (1, 8, 0):
        logger.info(f"{env_name}: pharokka {pharokka_version} already OK, skipping pharokka/mmseqs2/dnaapler fix.")
        return

    logger.info(f"{env_name}: pharokka {pharokka_version or 'unknown'} predates 1.8.0, applying fix.")
    run_live(["conda", "remove", "-n", env_name, "pharokka", "--force", "-y"], f"{env_name}: removing old conda pharokka")
    run_live(["conda", "install", "-n", env_name, "-y", "mmseqs2>=14"], f"{env_name}: upgrading mmseqs2")
    run_live(["conda", "run", "-n", env_name, "pip", "install", "pharokka==1.8.0"], f"{env_name}: installing pharokka 1.8.0")
    run_live(["conda", "run", "-n", env_name, "pip", "install", "dnaapler>=1.0.1"], f"{env_name}: upgrading dnaapler")

# PT_GNU_STACK program header type; PF_X is its executable-permission flag bit.
PT_GNU_STACK = 0x6474E551
PF_X = 0x1

# Clears the executable bit on an ELF file's GNU_STACK segment in place.
# Returns True if a change was made, False if it was already clear.
# (conda-forge's patchelf build here has no --clear-execstack option, and the
# `execstack` tool needs `prelink`, which isn't reliably installable either, so
# this patches the ELF program header directly instead of shelling out.)
def clear_execstack(path):
    with open(path, "r+b") as f:
        ident = f.read(16)
        if ident[:4] != b"\x7fELF" or ident[4] != 2 or ident[5] != 1:
            raise ValueError(f"{path} is not a little-endian 64-bit ELF file")

        f.seek(0)
        header = f.read(64)
        e_phoff, = struct.unpack_from("<Q", header, 32)
        e_phentsize, e_phnum = struct.unpack_from("<HH", header, 54)

        changed = False
        for i in range(e_phnum):
            offset = e_phoff + i * e_phentsize
            f.seek(offset)
            p_type, p_flags = struct.unpack("<II", f.read(8))
            if p_type == PT_GNU_STACK and (p_flags & PF_X):
                f.seek(offset + 4)
                f.write(struct.pack("<I", p_flags & ~PF_X))
                changed = True
        return changed

# pytorch 2.3.1's libtorch_cpu.so has the execstack bit set, which the kernel's
# security policy blocks. Clearing it in place fixes this (a version downgrade
# isn't viable here since medaka needs this torch series).
def fix_qppl_env2():
    env_name = "QPPL_env2"
    if not conda_env_exists(env_name):
        logger.warning(f"Conda environment '{env_name}' does not exist, skipping its fixes.")
        return

    returncode, torch_lib_path, _ = run_quiet([
        "conda", "run", "-n", env_name, "python", "-c",
        "import torch, os; print(os.path.join(os.path.dirname(torch.__file__), 'lib', 'libtorch_cpu.so'))"
    ])
    if returncode != 0 or not torch_lib_path or not os.path.isfile(torch_lib_path):
        logger.warning(f"{env_name}: could not locate libtorch_cpu.so, skipping execstack fix.")
        return

    try:
        if clear_execstack(torch_lib_path):
            logger.info(f"{env_name}: cleared execstack bit on {torch_lib_path}.")
        else:
            logger.info(f"{env_name}: execstack bit already clear on {torch_lib_path}.")
    except (OSError, ValueError) as e:
        logger.error(f"{env_name}: failed to clear execstack bit on {torch_lib_path}: {e}")

# phabox2 needs numpy<2 (yml pins 2.1.3) and pytorch 2.4.1 hits the same execstack
# crash as above; here it's fixed by downgrading torch instead, since phabox2
# doesn't hard-pin a torch version the way medaka does.
def fix_qppl_env3():
    env_name = "QPPL_env3"
    if not conda_env_exists(env_name):
        logger.warning(f"Conda environment '{env_name}' does not exist, skipping its fixes.")
        return

    numpy_version = _pip_version(env_name, "numpy")
    if not numpy_version or _version_tuple(numpy_version) >= (2, 0, 0):
        logger.info(f"{env_name}: pinning numpy<2 (currently {numpy_version or 'unknown'}).")
        run_live(["conda", "run", "-n", env_name, "pip", "install", "numpy<2"], f"{env_name}: installing numpy<2")
    else:
        logger.info(f"{env_name}: numpy {numpy_version} already OK.")

    torch_version = _pip_version(env_name, "torch")
    if torch_version != "2.1.2+cpu":
        logger.info(f"{env_name}: torch {torch_version or 'unknown'} needs the execstack-safe 2.1.2+cpu build.")
        run_live([
            "conda", "run", "-n", env_name, "pip", "install", "torch==2.1.2+cpu",
            "--index-url", "https://download.pytorch.org/whl/cpu"
        ], f"{env_name}: installing torch==2.1.2+cpu")
    else:
        logger.info(f"{env_name}: torch {torch_version} already OK.")

# Function to clone the vHULK repository if it isn't already present
def clone_vhulk(vhulk_dir):
    if os.path.isdir(os.path.join(vhulk_dir, ".git")):
        logger.info(f"vHULK repository already present at {vhulk_dir}, skipping clone.")
        return True

    os.makedirs(os.path.dirname(vhulk_dir.rstrip("/")) or "/", exist_ok=True)
    logger.info(f"Cloning vHULK into {vhulk_dir}...")
    return run_live(["git", "clone", VHULK_REPO_URL, vhulk_dir], "git clone vHULK")

# Main function to create/fix all conda environments and clone vHULK
def run_install_envs(vhulk_dir):
    log_path = os.path.join(REPO_ROOT, "install_envs.log")
    setup_logging(log_path)

    total_start_time = time.time()
    logger.info("=" * 60)
    logger.info("Installing/fixing QPPL conda environments")

    for env_name, yml_relpath in ENV_SPECS:
        create_env(env_name, yml_relpath)

    fix_qppl_env()
    fix_qppl_env2()
    fix_qppl_env3()

    clone_vhulk(vhulk_dir)

    total_elapsed = time.time() - total_start_time
    logger.info(f"Finished installing/fixing conda environments (elapsed time: {total_elapsed:.2f} seconds)")
    logger.info(f"vHULK location: {vhulk_dir} (set 'vhulk_location' in qppl.conf to this path if it differs from the default)")
