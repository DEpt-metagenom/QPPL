<p align="center">
  <img src="images/QPPL.png">
</p>

# QPPL (Quick Phage PipeLine)

This tool is designed to aid the work with phages, that were sequenced using Nanopore technology. QPPL is using multiple tools to create and characterize phage genome assemblies while every data stays on the user's computer.

## Fixes implemented by Tinhornfish

Goldrush/goldpolish deadlock — Fixed by writing a watchdog which detects zombie child processes then kills the entire process group with `SIGKILL`.

Added `.{basename}.done` cache files for reruns

Missing `-p` on mkdir for reruns

`gzip noadapter.fastq` silently failed because `.fastq.gz` already existed — ixed by detecting whether the noadapter file already exists and skipping porechop entirely if so

Missing `-f` flag on pharokka for reruns

`mv ~/vHULK/{basename}` into an existing directory caused outputs to land at `out_vhulk/{basename}/{basename}/predictions/`. Fixed by doing `rm -rf {dest}` before `mv`

### Conda Environment Fixes — MANUAL

Pharokka v1.7.3 → v1.8.0 — Old conda package pinned MMseqs2 to v13; pharokka v1.8.0 requires v14+. Fixed by `conda remove pharokka --force`, then `conda install mmseqs2>=14`, then `pip install pharokka==1.8.0`

dnaapler version — Pharokka v1.8.0 requires dnaapler ≥ v1.0.1 but v0.8.1 was installed. Fixed with `pip install "dnaapler>=1.0.1"`

NumPy 1.x/2.x incompatibility — phabox2 broke with numpy 2.x. Fixed with `pip install "numpy<2"` in both QPPL_env and QPPL_env3

PyTorch execstack crash — PyTorch 2.4.1 had the execstack bit set, which the kernel security policy blocked. Fixed by downgrading to `torch==2.1.2+cpu`

taxmyPHAGE database filename **caveman fix** — taxmyPHAGE v0.3.6 expected ICTV_2024.msh but the database file was named ICTV.msh. Fixed with `mv ICTV.msh ICTV_2024.msh`

## Usage

1. Clone the repository:

   `git clone https://github.com/DEpt-metagenom/QPPL.git`

2. Create required Conda environments from `envs/`:

   `conda env create -f envs/QPPL_env.yml`

   `conda env create -f envs/QPPL_env2.yml`

   `conda env create -f envs/QPPL_env3.yml`

   `conda env create -f envs/QPPL_env4.yml`

   `conda env create -f envs/vHULK_env.yml`

3. Clone the vHULK repository and install its database.

4. Install all databases to a specific location (CheckV, Pharokka, taxmyPHAGE, PhaBOX, vHULK).

5. Generate a default config:

   `python QPPL.py -gc`

6. Edit `qppl.conf` paths/tasks/parameters as needed.

7. Put input reads into `input_dir` (configured in `[General]`).

   Supported extensions: `.fastq`, `.fastq.gz`, `.fq`, `.fq.gz`.

8. Run:

   `python QPPL.py -c qppl.conf`

## Command-line options

`python QPPL.py -h`

Available options:

- `-h`, `--help`: show help menu
- `-v`, `--version`: show logo/version
- `-gc`, `--generate-config`: create default `qppl.conf`
- `-c`, `--config`: path to config file (default: `qppl.conf`)

## Tasks

`[General] task` can be one of:

- `readquality`: run read quality checks only
- `filter`: run filtering only
- `assembly`: run assembly workflow only
- `characterization`: run characterization workflow only
- `all`: run `filter` -> `assembly` -> `characterization`

Notes:

- `readquality` is standalone (not part of `all`).
- `assembly` expects filtered reads under `output_dir/filter_task`
- `characterization` expects assembly outputs under `output_dir/assembly_task`

## Configuration

Generate config with `python QPPL.py -gc`. Current structure:

- `[General]`
  - `input_dir`
  - `output_dir`
  - `task`

- `[Filter]`
  - `cut_adapter` (`y`/`n`)
  - `filter_quality`
  - `filter_length`
  - `porechop_threads`

- `[Assembly]`
  - `flye_mode`, `flye_meta`, `flye_threads`
  - `goldrush_genome_size`, `goldrush_min_read_length`, `goldrush_threads`
  - `raven_threads`
  - `shasta_read_length`, `shasta_threads`
  - `wtdbg2_genome_size`, `wtdbg2_read_length`, `wtdbg2_kmer_fsize`, `wtdbg2_kmer_psize`, `wtdbg2_min_similarity`, `wtdbg2_realign`, `wtdbg2_threads`
  - `polish_threads`
  - `checkv_db`, `checkv_threads`

- `[Characterization]`
  - `pharokka_threads`, `pharokka_db`
  - `taxmyphage_threads`, `taxmyphage_db`
  - `phabox_db`, `phabox_threads`
  - `vhulk_location`, `vhulk_threads`

## External dependencies and databases

QPPL runs tools inside dedicated Conda environments and expects the following resources to be installed/configured:

- CheckV database (`checkv_db`)
- Pharokka database (`pharokka_db`)
- taxmyPHAGE database (`taxmyphage_db`)
- PhaBOX database (`phabox_db`)
- vHULK installation and database (`vhulk_location`)

Also ensure `conda` is available in your shell PATH.


## How to Run

1. **[FIRST RUN ONLY]** Clone the QPPL repository and move into it.
2. **[FIRST RUN ONLY]** Create all Conda envs from `envs/`.
3. **[FIRST RUN ONLY]** Install databases for CheckV, Pharokka, taxmyPHAGE, PhaBOX.
4. **[FIRST RUN ONLY]** Clone the vHULK repository, move into it and install its database.
5. Run `python QPPL.py -gc` to generate the default config file.
6. Edit `qppl.conf` to reflect system resources, database locations, preferred parameters, selected task.
7. Place raw reads into `input_dir`, defined in `qppl.conf`.
8. Run `python QPPL.py` (uses `qppl.conf` by default), or run `python QPPL.py -c [somename].conf` to use a different config file.

## Directory Structure

Inside `output_dir`, QPPL writes:

- `readquality_task/` (if `task=readquality`)
- `filter_task/`
- `assembly_task/`
- `characterization_task/`
- `results/`

Main logs are written to:

- `readquality.log`
- `filter.log`
- `assembly.log`
- `characterization.log`

Final per-sample outputs are written to:

- `results/genomes/<sample>/` (selected genomes)
- `results/<sample>_final_results.csv`

## Citation

If you use this piece of software, please refer to the following presentation:
Z. Halász et al. (2025). An automated pipeline to genotype bacteriophages using third generation sequencing. ESCMID Global, ePoster E0539, Vienna
