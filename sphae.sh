#!/bin/bash

#SBATCH --job-name=sphae
#SBATCH --output=%x-%j.out.txt
#SBATCH --error=%x-%j.err.txt
#SBATCH --time=8:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=200G

#sphae install --db_dir sphae_databases
export VVDB=/users/bnalaga1/scratch/reference_db/sphae_databases/Pfam35.0/Pfam-A.hmm.gz
export CHECKVDB=/users/bnalaga1/scratch/reference_db/sphae_databases/checkv-db-v1.5
export PHAROKKADB=/users/bnalaga1/scratch/reference_db/sphae_databases/pharokka_db
export PHYNTENYDB=/users/bnalaga1/scratch/reference_db/sphae_databases/models
export PHOLDDB=/users/bnalaga1/scratch/reference_db/sphae_databases/phold

sphae run --input tests/data/illumina-subset --threads 32 -k
sphae run --input tests/data/nanopore-subset --sequencing longread --threads 64 -k
sphae run --input tests/data/nanopore-subset --sequencing longread --threads 64 -k --no_medaka  
sphae annotate --genome /users/bnalaga1/scratch/soil_metagenomes/Comb_phage_genomes --threads 32 -k --profile slurm
#singularity exec -B /home/nala0006/scratch/sphae/sphae/workflow/databases:/databases,/home/nala0006/scratch/sphae/tests/data/illumina-subset:/input,/home/nala0006/scratch/sphae/output:/output /home/nala0006/scratch/docker/sphae_v1.4.8.sif  sphae run --input /input --output /output

