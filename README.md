[![Edwards Lab](https://img.shields.io/badge/Bioinformatics-EdwardsLab-03A9F4)](https://edwards.flinders.edu.au)
[![DOI](https://zenodo.org/badge/403889262.svg)](https://zenodo.org/doi/10.5281/zenodo.8365088)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
![GitHub language count](https://img.shields.io/github/languages/count/linsalrob/spae)
[![](https://img.shields.io/static/v1?label=CLI&message=Snaketool&color=blueviolet)](https://github.com/beardymcjohnface/Snaketool)
![GitHub last commit (branch)](https://img.shields.io/github/last-commit/linsalrob/spae/main)
[![CI](https://github.com/linsalrob/spae/actions/workflows/testing.yml/badge.svg)](https://github.com/linsalrob/spae/actions/workflows/testing.yml)

# Sphae 
## Phage toolkit to detect phage candidates for phage therapy
<p align="center">
  <img src="sphae.png#gh-light-mode-only" width="300">
  <img src="sphaedark.png#gh-dark-mode-only" width="300">
</p>



**Overview**

This snakemake workflow was built using Snaketool [https://doi.org/10.1371/journal.pcbi.1010705], to assemble and annotate phage sequences. Currently, this tool is being developed for phage genomes. The steps include,

- Quality control that removes adaptor sequences, low-quality reads and host contamination (optional). 
- Assembly
- Contig quality checks; read coverage, viral or not, completeness, and assembly graph components. 
- Phage genome annotation'
- Annotation of the phage genome 
  
A complete list of programs used for each step is mentioned in the sphae.CITATION file. 

### Install 

**Pre-requisites**   
  - gcc
  - conda 
  - libgl1-mesa-dev (ubuntu- for Bandage)
  - libxcb-xinerama0 (ubuntu- for Bandage)

**Install**

Setting up a new conda environment 

```bash
conda create -n sphae python=3.11
conda activate sphae
conda install -n base -c conda-forge mamba #if you don't already have mamba installed
```

Steps for installing sphae workflow 

```bash
#clone sphae repository
git clone https://github.com/linsalrob/sphae.git

#move to sphae folder
cd sphae

#install sphae
pip install -e .

#confirm the workflow is installed by running the below command 
sphae --help
```

## Installing databases
Run command,

```
sphae install
```

  Install the databases to a directory, `sphae/workflow/databases`

  This workflow requires the 
  - Pfam35.0 database to run viral_verify for contig classification. 
  - CheckV database to test for phage completeness
  - Pharokka databases 
  - Phynteny models

This step takes approximately 1hr 30min to install and requires 9G of storage

## Running the workflow

The command `sphae run` will run QC, assembly and annotation

**Commands to run**

Only one command needs to be submitted to run all the above steps: QC, assembly and assembly stats

```bash
#For illumina reads, place the reads both forward and reverse reads to one directory
sphae run --input tests/data/illumina-subset --output example -k 

#For nanopore reads, place the reads, one file per sample in a directory
sphae run --input tests/data/nanopore-subset --sequencing longread --output example -k

#To run either of the commands on the cluster, add --profile slurm to the command. For instance here is the command for longreads/nanopore reads 
#Before running this below command, make sure have slurm config files setup, here is a tutorial, https://fame.flinders.edu.au/blog/2021/08/02/snakemake-profiles-updated 
sphae run --input tests/data/nanopore-subset --preprocess longread --output example --profile slurm -k
```

**Output**

Output is saved to `example/RESULTS` directory. In this directory, there will be four files 
  - Genome annotations in GenBank format (Phynteny output)
  - Genome in fasta format (either the reoriented to terminase output from Pharokka, or assembled viral contigs)
  - Circular visualization in `png` format (Pharokka output)
  - Genome summary file

Genome summary file includes the following information to help, 
  - Sample name
  - Length of the genome 
  - Coding density
  - If the assembled contig is circular or not (From the assembly graph)
  - Completeness (calculated from CheckV)
  - Contamination (calculated from CheckV)
  - Taxonomy accession ID (Pharokka output, searches the genome against INPHARED database using mash)
  - Taxa mash includes the number of matching hashes of the assembled genome to the accession ID/Taxa name. Higher the matching hash- more likely the genome is related to the taxa predicted
  - Gene searches:
    - Whether integrase is found (search for integrase gene in annotations)
    - Whether anti-microbial genes were found (Pharokka search against AMR database)
    - Whether any virulence factors were found (Pharokka search against virulence gene database)
    - Whether any CRISPR spacers were found (Pharokka search against MinCED database) 
 
## Issues and Questions

This is still a work in progress, so if you come across any issues or errors, report them under Issues. 

