[![Edwards Lab](https://img.shields.io/badge/Bioinformatics-EdwardsLab-03A9F4)](https://edwards.flinders.edu.au)
[![DOI](https://zenodo.org/badge/403889262.svg)](https://zenodo.org/doi/10.5281/zenodo.8365088)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

![GitHub language count](https://img.shields.io/github/languages/count/linsalrob/spae)
[![](https://img.shields.io/static/v1?label=CLI&message=Snaketool&color=blueviolet)](https://github.com/beardymcjohnface/Snaketool)
![GitHub last commit (branch)](https://img.shields.io/github/last-commit/linsalrob/spae/main)
[![CI](https://github.com/linsalrob/spae/actions/workflows/testing.yml/badge.svg)](https://github.com/linsalrob/spae/actions/workflows/testing.yml)

[![install with pip](https://img.shields.io/static/v1?label=Install%20with&message=PIP&color=success)](https://pypi.org/project/sphae/)
[![Pip Downloads](https://static.pepy.tech/badge/sphae)](https://www.pepy.tech/projects/sphae)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/sphae/README.html)
[![Bioconda Downloads](https://img.shields.io/conda/dn/bioconda/sphae)](https://img.shields.io/conda/dn/bioconda/sphae)
![Docker Pulls](https://img.shields.io/docker/pulls/npbhavya/sphae.svg)


# Sphae 
## Phage toolkit to detect phage candidates for phage therapy
<p align="center">
  <img src="img/sphae.png#gh-light-mode-only" width="300">
  <img src="img/sphaedark.png#gh-dark-mode-only" width="300">
</p>


**Overview**

The steps that sphae takes are shown here:
<p align="center">
  <img src="img/sphae_steps.png#gh-light-mode-only" width="300">
</p>

This snakemake workflow was built using Snaketool [https://doi.org/10.1371/journal.pcbi.1010705], to assemble and annotate phage sequences. Currently, this tool is being developed for phage genomes. The steps include,

- Quality control that removes adaptor sequences, low-quality reads and host contamination (optional). 
- Assembly
- Contig quality checks; read coverage, viral or not, completeness, and assembly graph components. 
- Phage genome annotation

**Cite Sphae: https://doi.org/10.1093/bioadv/vbaf004**

**If you are new to bioinformatics or running command line tools, here is a great tutorial to follow: https://github.com/AnitaTarasenko/sphae/wiki/Sphae-tutorial**

### Install 

**Pip install**

```bash
#creating a new envrionment
conda create -y -n sphae python=3.12
conda activate sphae
#install sphae 
pip install sphae
```

**Conda install** 

Setting up a new conda environment 

```bash
conda create -n sphae python
conda activate sphae
#if you don't already have mamba installed
conda install -n base -c conda-forge mamba
```

**Container Install**

We have two containers available, 
1. [Sphae v1.4.6 with databases](https://hub.docker.com/repository/docker/npbhavya/sphae)
   This is very large container, about 20.31 GB, so it may take a while to download and install.

   Here are the commands to download sphae container with databases
    ```
    TMPDIR=<where your tmpdir lives>
    IMAGEDIR-<where you want the image to live>
    
    singularity pull --tmpdir=$TMPDIR --dir $IMAGEDIR docker://npbhavya/sphae:latest
    singularity exec sphae_latest.sif sphae --help
    singularity exec sphae_latest.sif sphae run --help
    singularity exec sphae_latest.sif sphae install --help

    singularity exec -B <path/to/inputfiles>:/input,<path/to/output>:/output sphae_latest.sif sphae run --input /input --output /output
    ```
    
2. [Sphae v1.4.5 **without** databases](https://hub.docker.com/repository/docker/npbhavya/sphae)
   This version of sphae container does not include the databases, so they would have to be downloaded separately. The advantage of this is the container is smaller, so quick to donwnload and the databases can be downloaded separately. 

   You will still need to install the databases with `sphae install` as outlined below.

   ```
   TMPDIR=<where your tmpdir lives>
    IMAGEDIR-<where you want the image to live>
    
    singularity pull --tmpdir=$TMPDIR --dir $IMAGEDIR docker://npbhavya/sphae:latest
    #test if sphae is installed 
    singularity exec sphae_latest.sif sphae --help
    singularity exec sphae_latest.sif sphae run --help
    #mount the databases and input files to the image and run with a dataset
    singularity exec -B </path/to/databases>:/databases, <path/to/inputfiles>:/input,<path/to/output>:/output sphae_latest.sif sphae run --input /input --db_dir /databases --output /output
   ```
   
**Source install**

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
Run the below command,

```bash
#Installs the database to default directory, `sphae/workflow/databases`
sphae install

#Install database to specific directory
sphae install --db_dir <directory> 
```

  Install the databases to a directory, `sphae/workflow/databases`

  This workflow requires the 
  - Pfam35.0 database to run viral_verify for contig classification. 
  - CheckV database to test for phage completeness
  - Pharokka databases 
  - Phynteny models
  - Phold databases

This step requires ~17G of storage

If these databases are already installed, skip thsi step and instead set the envrionment variables pointing to the where these databases are installed

```bash
export VVDB=sphae/workflow/databases/Pfam35.0/Pfam-A.hmm.gz
export CHECKVDB=sphae/workflow/databases/checkv-db-v1.5
export PHAROKKADB=sphae/workflow/databases/pharokka_db
export PHYNTENYDB=sphae/workflow/databases/phynteny_models_zenodo
export PHOLDDB=sphae/workflow/databases/phold
```

## Running the workflow

Sphae is developed to be modular: 
- `sphae run` will run QC, assembly and annotation
- `sphae annotate` will run only annotation steps
  
**Commands to run**

Only one command needs to be submitted to run all the above steps: QC, assembly and assembly stats

```bash
#For illumina reads, place the reads both forward and reverse reads to one directory
#Make sure the fastq reads are saved as {sample_name}_R1.fastq and {sample_name}_R2.fastq or with extensions {sample_name}_R1.fastq.gz
sphae run --input tests/data/illumina-subset --output example -k --use-conda --conda-frontend mamba

#For nanopore reads, place the reads, one file per sample in a directory
sphae run --input tests/data/nanopore-subset --sequencing longread --output example -k --use-conda --conda-frontend mamba

#For newer ONT sequencing data where polishing is not required, run the command
sphae run --input tests/data/nanopore-subset --sequencing longread --output example -k --no_medaka --use-conda --conda-frontend mamba

#To run either of the commands on the cluster, add --executor slurm to the command. There is a little bit of setup to do here.
#Setup a ~/.config/snakemake/slurm/config.yaml file - https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/slurm.html#advanced-resource-specifications
#I may have set this workflow to run only slurm right now, will make it more generic soon.
sphae run --input tests/data/nanopore-subset --preprocess longread --output example --profile slurm -k --threads 16 --use-conda --conda-frontend mamba

```

**Command to run only annotation steps and phylogenetic trees**
This step reruns 
   - Pharokka, Phold, Phynteny
   - Phylogenetic tree with terminase large subunit, portal protein
   
```bash
#the genomes directory has the already assembled complete genomes
sphae annotate --genome <genomes directory> --output example -k --use-conda --conda-frontend mamba
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
  - If the assembled contig is circular or not (From checkv)
  - Completeness (calculated from CheckV)
  - Contamination (calculated from CheckV)
  - Taxonomy accession ID (Pharokka output, searches the genome against INPHARED database using mash)
  - Taxa mash includes the number of matching hashes of the assembled genome to the accession ID/Taxa name. Higher the matching hash- more likely the genome is related to the taxa predicted
  - Gene searches:
    - Whether integrase is found (search for integrase gene in annotations)
    - Whether anti-microbial genes were found (Phold and Pharokka search against AMR database)
    - Whether any virulence factors were found (Pharokka search against virulence gene database)
    - Whether any CRISPR spacers were found (Pharokka search against MinCED database) 

### FAQ
1. **"Failed during assembly":**
   - This message indicates that the assembly process was unsuccessful. It suggests that the assembler could not generate contigs, which are contiguous sequences of DNA, typically representing segments of a genome. 
   - To confirm this, you can check the logs located at `sphae.out/PROCESSING/assembly/flye/<sample name>/assembly_info.txt` or `sphae.out/PROCESSING/assembly/megahit/<sample name>/log`. These logs should provide details about the error or the step at which the assembly failed.
   - One possible reason for this failure could be insufficient genome coverage, meaning that there were not enough sequencing reads to accurately assemble the genome.

2. **"Genome includes multiple contigs, fragmented":**
   - This message indicates that the assembly generated numerous short fragments (contigs) instead of a single, contiguous sequence representing a nearly complete phage genome. 
   - You can verify this by examining the file `sphae.out/PROCESSING/assembly/flye/<sample name>-assembly-stats_flye.csv` or `sphae.out/PROCESSING/assembly/megahit/<sample name>-assembly-stats_megahit.csv`.
   - Each row in these tables represents a contig along with its characteristics. If none of the contigs are identified as viral and do not meet a certain completeness threshold (e.g., greater than 70% completeness), it suggests that the assembly consists of fragmented contigs.
   - Fragmented contigs make it challenging to accurately identify genes. To address this issue, you may need to resequence the phages for better coverage or try using different assembly algorithms.

3. **"Good genome coverage but still encountering assembly issues":**
   - If you have adequate genome coverage but still face assembly problems, you may consider adjusting the subsampling step in sphae. This step involves randomly selecting a subset of reads to reduce the computational burden.
   - To modify the subsampling parameters, navigate to the `config/config.yaml` file and update the line under `subsample` section, for example:
     ```
     subsample:
         --bases 1000M
     ```
   - Increase or decrease the number of bases (e.g., `1000M` for 1000 megabases) based on your requirements.
   - After making the changes, rerun sphae and ensure that the updated subsampling parameters are reflected in the `sphae.out/sphae.config.yaml` file.

4. **"What does 'No integrases found ...but Phynteny predicted a few unknown function genes to have some similarity with integrase genes but with low confidence. Maybe a false positive or a novel integrase gene' mean?"**
   This message indicates that while no integrase genes were explicitly identified, the analysis detected certain genes that exhibited similarities to integrase genes. However, these genes were associated with low confidence scores, suggesting a possibility of being false positives or potentially representing novel integrase genes.
   
   [Phynteny](https://github.com/susiegriggo/Phynteny), the tool used for this prediction, assigns a confidence score to each gene prediction. If this score falls below a certain threshold (typically 90%), the gene remains classified as having an unknown function. To further investigate these genes, advanced techniques such as folding using tools like [ColabFold](https://github.com/sokrypton/ColabFold) and [Foldseek](https://github.com/steineggerlab/foldseek) can be employed. Analyzing the structure of these genes may provide additional insights into their functionality and potential role in biological processes.

5. **How do I visualize the phages and gene annotations?**
   To visualize the phages and gene annotations, I recommend using [Clinker](https://github.com/gamcil/clinker). First, gather all the sample genbank files from `sphae.out/RESULTS` and place them in a new directory. Then, execute the clinker command to generate clinker plots, which compare the genes in each genome to each other.
   
   Additionally, for enhanced visualization, consider running [dnaapler](https://github.com/gbouras13/dnaapler) on the genomes in fasta format obtained from 
   `sphae.out/RESULTS`. This step generates reoriented phages that start with terminase genes. Pharokka -> Phold -> Phynteny has to be rerun, and the resulting genbank files can be used for visualization. To perform the annotation steps, run the command 
   `sphae annotate --input <reoriented genomes from dnaapler in fasta format directory>`
   
   Please note that dnaapler may fail if terminase genes are not found, particularly when working with novel phages. The reason these steps haven't been added to sphae. If you encounter any challenges during this process, please feel free to leave an issue, and I'll provide improved documentation to assist you further with the command on how to install and run the command different commands. 

6. **Where are the intermediate files being saved?**
   These files are being saved in `sphae.out/PROCESSING`. If you need more information on the file structure here, or have ideas of better organization then leave an issue and I will make a note to have more documentation. 

7. **Just run annotation on already assembled genomes?**
   
    `sphae annotate --input <input genomes>`
    This command runs only Pharokka, Phold and Phynteny to annotate the assembled genomes. The results are saved to a new directory labeled `sphae.out/annotation`. 

    Note: Currently, Sphae runs Phold in CPU mode, but efforts are underway to support Phold GPU mode for faster processing of this step.

8. How to change the number of base pairs to subsample for a sample?
    Run the command `sphae config`
    This copies the config file within the workflow to the current directory. Open this file and update the line `bases: 10000000` to for instance `bases: 300000`
    Then run sphae run with the command `sphae run --input tests/data/illumina-subset --output example -k --config <path to the config file with the change>`
   
   
## Citation
To cite sphae, doi: https://doi.org/10.1101/2024.11.18.624194

## Issues and Questions

If you come across any issues or errors, report them under [Issues](https://github.com/linsalrob/sphae/issues). 


