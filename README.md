<p align="center">
  <img src="images/QPPL.png">
</p>

# QPPL (Quick Phage PipeLine)

This tool is designed to aid the work with phages, that were sequenced using Nanopore technology. QPPL is using multiple tools to create and characterize phage genome assemblies while every data stays on the user's computer.

## Usage

- Clone this repository.

- Install all conda environments from `envs/` using `conda env create -f [env.yml]`.

- Clone the vHULK repository and install its database.

- Install all databases to a specific location (CheckV, Pharokka, PhaBOX, vHULK).

- Generate a config file using `python QPPL.py -gc` or `python QPPL.py --generate-config`.

  - Task by default is set to `all`, so it will run all automatized tasks: `filter`, `assembly`, `characterization`. These tasks can also be run separately.
 
  - The `readquality` task is not part of the automatized pipeline, but can be run by changing the task. It is useful for manual read quality assessment, this way the config file's filter part can be customized to fit the reads more.

  - Change the CheckV, Pharokka and PhaBOX databases listed in the config file to the correct location.

  - Change vHULK location in the config file to the correct location. 

- Create a directory containing the raw reads. Format can be `.fastq/.fq` or `.fastq.gz/.fq.gz`.

- Run QPPL using `python QPPL.py -c qppl.conf`.


## Quick Run

- (FIRST RUN ONLY) Clone the QPPL repository and move into it.

- (FIRST RUN ONLY) Install all conda environments.

- Generate a config file.

- (FIRST RUN ONLY) Install all databases where the config file says.

- (FIRST RUN ONLY) Clone the vHULK repository where the config file says, move into it and install its database.

- Make sure that the system can handle the thread values set by default. If not, please lower them. 

- Create an input directory called `input` in the QPPL directory and put all the raw data into it.

- Run QPPL.

## Directory Structure

-  The `output_dir` contains the following directories:

  - `readquality_task`

  - `filter_task`

  - `assembly_task`

  - `characterization_task`

  - `results`

- The `results` directory contains the final results: all phage genomes in the `genomes` subdirectory, and the characterization results as a `.csv`.

## Citation

If used, please refer to this Github repository.
