[![Edwards Lab](https://img.shields.io/badge/Bioinformatics-EdwardsLab-03A9F4)](https://edwards.flinders.edu.au)
[![DOI](https://zenodo.org/badge/403889262.svg)](https://zenodo.org/doi/10.5281/zenodo.8365088)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
![GitHub language count](https://img.shields.io/github/languages/count/linsalrob/spae)
[![](https://img.shields.io/static/v1?label=CLI&message=Snaketool&color=blueviolet)](https://github.com/beardymcjohnface/Snaketool)
![GitHub last commit (branch)](https://img.shields.io/github/last-commit/linsalrob/spae/main)
[![CI](https://github.com/linsalrob/spae/actions/workflows/testing.yml/badge.svg)](https://github.com/linsalrob/spae/actions/workflows/testing.yml)
[![Upload Python Package](https://github.com/linsalrob/spae/actions/workflows/python-publish.yml/badge.svg)](https://github.com/linsalrob/spae/actions/workflows/python-publish.yml)

# Sphae 
## Phage toolkit to detect phage candidates for phage therapy
<p align="center">
  <img src="sphae.png#gh-light-mode-only" width="300">
  <img src="sphaedark.png#gh-dark-mode-only" width="300">
</p>



**Overview**

This snakemake workflow was built using Snaketool [https://doi.org/10.1371/journal.pcbi.1010705], to assemble and annotate phage sequences. Currently this tool is being developed for phage genomes. The steps include,

- Quality control that removes adaptor sequences, low quality reads and host contimanination (optional). 
- Assembly
- Contig quality checks; read coverage, viral or not, completeness, and assembly graph components. 
- Phage genome annotation'
- Annotation of the phage genome 
  
Complete list of programs used for each step is mention in the sphae.CITATION file. 

### Install 

**Pre-requisites**   
  - gcc
  - conda 
  - libgl1-mesa-dev (ubuntu- for Bandage)
  - libxcb-xinerama0 (ubuntu- for Bandage)

**Install**
Setting up a new conda environment 

    conda create -n sphae python=3.11
    conda activate sphae
    conda install -n base -c conda-forge mamba #if you dont already have mamba installed

Steps for installing sphae workflow 

    git clone https://github.com/linsalrob/sphae.git
    cd sphae
    pip install -e .
    #confirm the workflow is installed by running the below command 
    sphae --help

## Installing databases
Run command,

    sphae install

  Install the databases to a directory, sphae/workflow/databases

  This workflow requires the 
  - Pfam35.0 database to run viral_verify for contig classification. 
  - CheckV database to test for phage completeness
  - Pharokka databases 
  - Phyteny models

This step takes approximately 1hr 30min to install, and requires 9G of storage

## Running the workflow
The command `sphae run` will run QC, assembly and annoation

**Commands to run**
Only one command needs to be submitted to run all the above steps: QC, assembly and assembly stats

    #For illumina reads, place the reads both forward and reverse reads to one directory
    sphae run --input tests/data/illumina-subset --output example

    #For nanopore reads, place the reads, one file per sample in a directory
    sphae run --input tests/data/nanopore-subset --sequencing longread --output example 

    #To run either of the commands on the cluster, add --profile slurm to the command. For instance here is the command for longreads/nanopore reads 
    #Before running this below command, makse sure have slurm config files setup, here is a tutorial, https://fame.flinders.edu.au/blog/2021/08/02/snakemake-profiles-updated 
    sphae run --input tests/data/nanopore-subset --preprocess longread --output example --profile slurm 

**Output**
- Assmbled phage genome saved to **"{outut-directory}/genome/{sample}/{sample}.fasta**
- Annotations of the phage genome are saved to **"{outut-directory}/pharokka/phynteny/phynteny.gbk"**
 
## Issues and Questions

This is still a work in progress, so if you come across any issues or errors, report them under Issues. 

