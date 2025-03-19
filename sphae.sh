#!/bin/bash

#SBATCH --job-name=sphae
#SBATCH --mail-type=ALL
#SBATCH --output=%x-%j.out.txt
#SBATCH --error=%x-%j.err.txt
#SBATCH --time=1-0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=200G
#SBATCH --partition=high-capacity
#SBATCH --qos=hc-concurrent-jobs

#sphae install
export VVDB=/home/nala0006/scratch/sphae/sphae/workflow/databases/Pfam35.0/Pfam-A.hmm.gz
export CHECKVDB=/home/nala0006/scratch/sphae/sphae/workflow/databases/checkv-db-v1.5
export PHAROKKADB=/home/nala0006/scratch/sphae/sphae/workflow/databases/pharokka_db
export PHYNTENYDB=/home/nala0006/scratch/sphae/sphae/workflow/databases/phynteny_models_zenodo
export PHOLDDB=/home/nala0006/scratch/sphae/sphae/workflow/databases/phold

#sphae run --input tests/data/illumina-subset --threads 32 -k
#sphae run --input tests/data/nanopore-subset --sequencing longread --threads 64 -k
#sphae run --input tests/data/nanopore-subset --sequencing longread --threads 64 -k --no_medaka  
sphae annotate --genome tests/data/genome --threads 64 -k
