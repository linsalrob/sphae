#!/bin/bash

#SBATCH --job-name=pharokka
#SBATCH --mail-type=ALL
#SBATCH --output=%x-%j.out.txt
#SBATCH --error=%x-%j.err.txt
#SBATCH --time=1-0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=500G
#SBATCH --partition=high-capacity
#SBATCH --qos=hc-concurrent-jobs

#create_custom_hmm.py -i /home/nala0006/scratch/spae/phage_genome_assembly/workflow/databases/pharokka_db/each_family_msa -o dbAPIS -p dbAPIS
pharokka.py -i ../spae-example/PROCESSING/genome/nanoreads-sr/nanoreads_genome.fasta -o ../spae-example/test-pharokka -d phage_genome_assembly/workflow/databases/pharokka_db -t 16 -f -p test-pharokka --dnaapler --custom_hmm phage_genome_assembly/workflow/databases/pharokka_db/dbAPIS/dbAPIS.h3m
#pharokka.py -i ../spae-example/PROCESSING/genome/nanoreads-sr/nanoreads_genome.fasta -o ../spae-example/PROCESSING/pharokka/nanoreads-sr -d /scratch/user/nala0006/spae/phage_genome_assembly/workflow/databases/pharokka_db -t 64 -f -p nanoreads --dnaapler 
