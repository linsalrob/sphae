#!/bin/bash

#SBATCH --job-name=sphae
#SBATCH --mail-type=ALL
#SBATCH --output=%x-%j.out.txt
#SBATCH --error=%x-%j.err.txt
#SBATCH --time=1-0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=500G
#SBATCH --partition=high-capacity
#SBATCH --qos=hc-concurrent-jobs

#sphae install
#sphae run --input tests/data/illumina-subset --threads 64 -k 
#sphae run --input /home/nala0006/scratch/Achromobacter/phage/fastq --sequencing longread --threads 64 -k --output /home/nala0006/scratch/Achromobacter/phage/sphae-outv1.4.4
#sphae run --input tests/data/nanopore-subset --sequencing longread --threads 64 -k --no_medaka 
sphae annotate --genome /home/nala0006/scratch/Achromobacter/phage/fasta2 --threads 64 --output /home/nala0006/scratch/Achromobacter/phage/2089-annot-v1.4.4

