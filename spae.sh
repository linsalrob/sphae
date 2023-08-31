#!/bin/bash
#SBATCH --job-name=spae
#SBATCH --output=%x-%j.out.txt
#SBATCH --error=%x-%j.err.txt
#SBATCH --time=1-0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=500G
#SBATCH --partition=high-capacity
#SBATCH --qos=hc-concurrent-jobs

spae run --input test/illumina-subset --output ../example-ill -k 
