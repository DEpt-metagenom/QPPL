<p align="center">
  <img src="images/QPPL.png">
</p>

# QPPL (Quick Phage PipeLine)

This tool is designed to aid the work with phages, that were sequenced using Nanopore technology. QPPL is using multiple tools to create and characterize phage genome assemblies while every data stays on the user's computer.

## Usage

1. Clone the repository:

   `git clone https://github.com/DEpt-metagenom/QPPL.git`

2. Create all required Conda environments, and clone vHULK:

   `python QPPL.py -ie`

   (Pass `--vhulk-dir <path>` to clone vHULK somewhere other than the default `~/vHULK`.)

3. Download the CheckV, Pharokka, taxmyPHAGE and PhaBOX databases into `Databases/`:

   `python QPPL.py -dd`

   (Pass `--databases-dir <path>` to use a different location.) Both `-ie` and `-dd`
   are safe to re-run — anything already installed/downloaded is skipped.

4. Generate a default config:

   `python QPPL.py -gc`

5. Edit `qppl.conf` paths/tasks/parameters as needed (point `checkv_db`, `pharokka_db`,
   `taxmyphage_db`, `phabox_db` at the paths under `Databases/`, and `vhulk_location`
   at the directory used in step 2).

6. Put input reads into `input_dir` (configured in `[General]`).

   Supported extensions: `.fastq`, `.fastq.gz`, `.fq`, `.fq.gz`.

7. Run:

   `python QPPL.py -c qppl.conf`

## Command-line options

`python QPPL.py -h`

Available options:

- `-h`, `--help`: show help menu
- `-v`, `--version`: show logo/version
- `-gc`, `--generate-config`: create default `qppl.conf`
- `-c`, `--config`: path to config file (default: `qppl.conf`)
- `-ie`, `--install-envs`: create (or fix) all conda environments needed by QPPL, and clone vHULK
- `--vhulk-dir`: destination directory to clone vHULK into (used with `-ie`, default `~/vHULK`)
- `-dd`, `--download-databases`: download the CheckV, Pharokka, taxmyPHAGE and PhaBOX databases into `Databases/`
- `--databases-dir`: destination directory for downloaded databases (used with `-dd`, default `Databases`)

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
  
`python QPPL.py -ie` and `python QPPL.py -dd` (see Usage above) **to set all of these up automatically**

Ensure `conda` is available in your shell PATH.

## How to Run

1. **[FIRST RUN ONLY]** Clone the QPPL repository and move into it.
2. **[FIRST RUN ONLY]** Run `python QPPL.py -ie` to create/fix all Conda envs and clone vHULK.
3. **[FIRST RUN ONLY]** Run `python QPPL.py -dd` to download the CheckV, Pharokka, taxmyPHAGE and PhaBOX databases into `Databases/`.
4. Run `python QPPL.py -gc` to generate the default config file.
5. Edit `qppl.conf` to reflect system resources, database locations, preferred parameters, selected task.
6. Place raw reads into `input_dir`, defined in `qppl.conf`.
7. Run `python QPPL.py` (uses `qppl.conf` by default), or run `python QPPL.py -c [filename].conf` to use a different config file.

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
