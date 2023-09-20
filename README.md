# Spae 
## Phage toolkit to detect phage candidates for phage therapy
<p align="center">
  <img src="spaefinal.png#gh-light-mode-only" width="400">
  <img src="spaedark.png#gh-dark-mode-only" width="400">
</p>
  
**Overview**

This snakemake workflow was built using Snaketool [https://doi.org/10.1371/journal.pcbi.1010705], to assemble and annotate phage sequences. Currently this tool is being developed for phage genomes. THe steps include 
- QC using Trimanammi [https://github.com/beardymcjohnface/Trimnami]
- Assembly, either SPAdes [https://github.com/ablab/spades] or Flye [https://github.com/fenderglass/Flye]
- Contig quality checks 
    - read coverage using Koverage [https://github.com/beardymcjohnface/Koverage]
    - verify they are viral contigs using ViralVerify [https://github.com/ablab/viralVerify]
    - Completeness using CheckV [https://bitbucket.org/berkeleylab/CheckV]
    - assembly graph compnents check - internal script 
  This results in determining if the phage genome was assembled 
- Annotation of the phage genome using Pharokka [https://github.com/gbouras13/pharokka]
    - To improve annotations, the results are run through Phynteny [https://github.com/susiegriggo/Phynteny] 

### Install 

**Pre-requisites**   
  - gcc
  - conda 
  - libgl1-mesa-dev (ubuntu- for Bandage)
  - libxcb-xinerama0 (ubuntu- for Bandage)

**Install**
Setting up a new conda environment 

    conda create -n spae python=3.11
    conda activate spae
    conda install -n base -c conda-forge mamba #if you dont already have mamba installed

Steps for installing spae workflow 

    git clone https://github.com/linsalrob/spae.git
    cd spae
    pip install -e .
    #confirm the workflow is installed by running the below command 
    spae --help

## Installing databases
Run command,

    spae install

  Install the databases to a directory, phage_genome_assembly/workflow/databases

  This workflow requires the 
  - Pfam35.0 database to run viral_verify for contig classification. 
  - CheckV database to test for phage completeness
  - Pharokka databases 
  - Phyteny models

This step takes approximately 1hr 30min to install

## Running the workflow
The command `spae run` will run QC, assembly and annoation

**Commands to run**
Only one command needs to be submitted to run all the above steps: QC, assembly and assembly stats

    #For illumina reads, place the reads both forward and reverse reads to one directory
    spae run --input test/illumina-subset --output example

    #For nanopore reads, place the reads, one file per sample in a directory
    spae run --input test/nanopore-subset --sequencing longread --output example 

    #To run either of the commands on the cluster, add --profile slurm to the command. For instance here is the command for longreads/nanopore reads 
    #Before running this below command, makse sure have slurm config files setup, here is a tutorial, https://fame.flinders.edu.au/blog/2021/08/02/snakemake-profiles-updated 
    spae run --input test/nanopore-subset --preprocess longread --output example --profile slurm 

**Output**
- Assmbled phage genome saved to **"{outut-directory}/genome/{sample}/{sample}.fasta**
- Annotations of the phage genome are saved to **"{outut-directory}/pharokka/phynteny/phynteny.gbk"**

**Issues and Questions**

This is still a work in progress, so if you come across any issues or errors, report them under Issues. 

