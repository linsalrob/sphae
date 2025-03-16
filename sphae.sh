#!/bin/bash

#SBATCH --job-name=sphae
#SBATCH --mail-type=ALL
#SBATCH --output=%x-%j.out.txt
#SBATCH --error=%x-%j.err.txt
#SBATCH --time=1-0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=200G
##SBATCH --partition=high-capacity
##SBATCH --qos=hc-concurrent-jobs

#sphae install
export VVDB=sphae/workflow/databases/Pfam35.0/Pfam-A.hmm.gz
export CHECKVDB=sphae/workflow/databases/checkv-db-v1.5
export PHAROKKADB=sphae/workflow/databases/pharokka_db
export PHYNTENYDB=sphae/workflow/databases/phynteny_models_zenodo
export PHOLDDB=sphae/workflow/databases/phold

#sphae run --input tests/data/illumina-subset --threads 32 -k --use-conda --conda-frontend mamba 
#sphae run --input tests/data/nanopore-subset --sequencing longread --threads 64 -k --use-conda --conda-frontend mamba 
sphae run --input tests/data/nanopore-subset --sequencing longread --threads 64 -k --no_medaka --use-conda --conda-frontend mamba 
sphae annotate --genome tests/data/genome --threads 64 -k --use-conda --conda-frontend mamba
