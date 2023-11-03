#!/bin/bash

#SBATCH --job-name=install-db
#SBATCH --mail-type=ALL
#SBATCH --output=%x-%j.out.txt
#SBATCH --error=%x-%j.err.txt
#SBATCH --time=1-0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=100G
#SBATCH --partition=high-capacity
#SBATCH --qos=hc-concurrent-jobs

time spae install
