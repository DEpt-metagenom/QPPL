import os
import glob
import shutil
import logging
import time
import urllib.request

from modules.setup import REPO_ROOT, conda_env_exists, run_live, run_quiet

# Call logger
logger = logging.getLogger("QPPL")

PHABOX_DB_URL = "https://github.com/KennthShang/PhaBOX/releases/download/v2/phabox_db_v2_2.zip"

VMR_URL = "https://ictv.global/vmr/current"
# ICTV's server 403s the default User-Agent used by both urllib and the `wget`
# pip package (which taxmyphage's own downloader relies on), but allows a
# browser-like one through to a 307 redirect -> 200. Verified with curl.
BROWSER_USER_AGENT = (
    "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 "
    "(KHTML, like Gecko) Chrome/124.0.0.0 Safari/537.36"
)

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

# Function to mark a database directory as fully downloaded
def _mark_done(dest):
    open(os.path.join(dest, ".done"), "w").close()

# Neither checkv's own downloader nor the PhaBOX zip has any built-in "already
# installed" check, so re-running would otherwise re-download from scratch every
# time. `expected_globs` are glob patterns (relative to dest) that only exist once
# a database is actually present; if any matches, treat the database as already
# downloaded and write the `.done` marker retroactively.
def _already_present(dest, expected_globs):
    if os.path.exists(os.path.join(dest, ".done")):
        return True
    if any(glob.glob(os.path.join(dest, pattern)) for pattern in expected_globs):
        _mark_done(dest)
        return True
    return False

# CheckV's own downloader extracts into a versioned subfolder (checkv-db-vX.Y/),
# but checkv_db/qppl.conf expects the flat directory, so flatten it afterwards.
def download_checkv(dest, env_name="QPPL_env"):
    if _already_present(dest, ["genome_db", "hmm_db"]):
        logger.info(f"CheckV database already present at {dest}, skipping.")
        return

    if not conda_env_exists(env_name):
        logger.error(f"Conda environment '{env_name}' does not exist. Run 'python QPPL.py -ie' first. Skipping CheckV database.")
        return

    os.makedirs(dest, exist_ok=True)
    logger.info(f"Downloading CheckV database into {dest}...")
    if not run_live(["conda", "run", "-n", env_name, "checkv", "download_database", dest], "checkv download_database"):
        return

    versioned_dirs = [d for d in glob.glob(os.path.join(dest, "checkv-db-v*")) if os.path.isdir(d)]
    for versioned_dir in versioned_dirs:
        for entry in os.listdir(versioned_dir):
            shutil.move(os.path.join(versioned_dir, entry), os.path.join(dest, entry))
        os.rmdir(versioned_dir)

    _mark_done(dest)
    logger.info(f"CheckV database ready at {dest}.")

# Pharokka's install_databases.py already knows the right zenodo URL/md5 for the
# installed pharokka version and self-checks whether the DB is already present.
def download_pharokka(dest, env_name="QPPL_env"):
    if _already_present(dest, ["pharokka_v*_databases"]):
        logger.info(f"Pharokka database already present at {dest}, skipping.")
        return

    if not conda_env_exists(env_name):
        logger.error(f"Conda environment '{env_name}' does not exist. Run 'python QPPL.py -ie' first. Skipping Pharokka database.")
        return

    os.makedirs(dest, exist_ok=True)
    logger.info(f"Downloading Pharokka database into {dest}...")
    if not run_live(["conda", "run", "-n", env_name, "install_databases.py", "-o", dest], "install_databases.py"):
        return

    _mark_done(dest)
    logger.info(f"Pharokka database ready at {dest}.")

# taxmyPHAGE's mash sketch is expected under both names depending on version/caller
# (ICTV.msh vs ICTV_2024.msh, see taxmyphage/__main__.py) — they're identical in
# content, just named differently, so keep both present via a copy either way.
def _ensure_ictv_alias(dest):
    old_msh = os.path.join(dest, "ICTV.msh")
    new_msh = os.path.join(dest, "ICTV_2024.msh")
    if os.path.exists(new_msh) and not os.path.exists(old_msh):
        logger.info(f"Copying {new_msh} to {old_msh}.")
        shutil.copy2(new_msh, old_msh)
    elif os.path.exists(old_msh) and not os.path.exists(new_msh):
        logger.info(f"Copying {old_msh} to {new_msh}.")
        shutil.copy2(old_msh, new_msh)

# taxmyphage's own VMR downloader (taxmyphage/download_check.py) uses the `wget`
# pip package with no custom headers, so ICTV's server always 403s it, install
# never gets further. Pre-fetch VMR.xlsx ourselves with a browser User-Agent so
# taxmyphage's check_VMR() finds it already present and skips its own downloader.
def _prefetch_vmr(dest):
    vmr_path = os.path.join(dest, "VMR.xlsx")
    if os.path.exists(vmr_path):
        return True

    logger.info(f"Pre-downloading ICTV VMR.xlsx into {dest} (taxmyphage's own downloader is blocked by ICTV's server)...")
    request = urllib.request.Request(VMR_URL, headers={"User-Agent": BROWSER_USER_AGENT})
    try:
        with urllib.request.urlopen(request, timeout=60) as response, open(vmr_path, "wb") as out_file:
            shutil.copyfileobj(response, out_file)
        return True
    except OSError as e:
        logger.error(f"Failed to pre-download VMR.xlsx from {VMR_URL}: {e}")
        if os.path.exists(vmr_path):
            os.remove(vmr_path)
        return False

def download_taxmyphage(dest, env_name="QPPL_env2"):
    if _already_present(dest, ["ICTV_2024.msh", "ICTV.msh"]):
        logger.info(f"taxmyPHAGE database already present at {dest}, skipping.")
        _ensure_ictv_alias(dest)
        return

    if not conda_env_exists(env_name):
        logger.error(f"Conda environment '{env_name}' does not exist. Run 'python QPPL.py -ie' first. Skipping taxmyPHAGE database.")
        return

    os.makedirs(dest, exist_ok=True)
    if not _prefetch_vmr(dest):
        return

    logger.info(f"Downloading taxmyPHAGE database into {dest}...")
    if not run_live(["conda", "run", "-n", env_name, "taxmyphage", "install", "--db_folder", dest], "taxmyphage install"):
        return

    _ensure_ictv_alias(dest)

    # taxmyphage calls bare `sys.exit()` on internal failures (e.g. a blocked
    # download), which exits 0 and makes run_live report success even though
    # nothing was installed. Verify the actual artifact exists before trusting that.
    if not glob.glob(os.path.join(dest, "ICTV_2024.msh")) and not glob.glob(os.path.join(dest, "ICTV.msh")):
        logger.error(f"taxmyphage install exited successfully but produced no mash index in {dest}; not marking as done.")
        return

    _mark_done(dest)
    logger.info(f"taxmyPHAGE database ready at {dest}.")

# PhaBOX has no bundled downloader; its database is a plain zip on GitHub releases.
def download_phabox(dest, url=PHABOX_DB_URL):
    if _already_present(dest, ["phabox_db_v2_2"]):
        logger.info(f"PhaBOX database already present at {dest}, skipping.")
        return

    os.makedirs(dest, exist_ok=True)
    zip_path = os.path.join(dest, "phabox_db_v2_2.zip")
    logger.info(f"Downloading PhaBOX database into {dest}...")
    if not run_live(["wget", "-c", url, "-O", zip_path], "wget phabox_db_v2_2.zip"):
        return

    if not run_live(["unzip", "-o", zip_path, "-d", dest], "unzip phabox_db_v2_2.zip"):
        return

    os.remove(zip_path)
    _mark_done(dest)
    logger.info(f"PhaBOX database ready at {dest}.")

# Main function to download all databases used by the characterization/assembly steps
def run_download_databases(databases_dir):
    setup_logging(os.path.join(REPO_ROOT, "databases.log"))

    total_start_time = time.time()
    logger.info("=" * 60)
    logger.info(f"Downloading QPPL databases into {databases_dir}")

    os.makedirs(databases_dir, exist_ok=True)

    download_checkv(os.path.join(databases_dir, "CheckV"))
    download_pharokka(os.path.join(databases_dir, "Pharokka"))
    download_taxmyphage(os.path.join(databases_dir, "taxmyPHAGE"))
    download_phabox(os.path.join(databases_dir, "PhaBOX"))

    total_elapsed = time.time() - total_start_time
    logger.info(f"Finished downloading databases (elapsed time: {total_elapsed:.2f} seconds)")
