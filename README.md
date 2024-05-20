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
- Phage genome annotation
  
A complete list of programs used for each step is mentioned in the `sphae.CITATION` file. 

### Install 

**Pip install**

```bash
pip install sphae
```

**Conda install** 
```bash
conda install sphae
```
**Source Install**

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
Run the below command,

```bash
#Installs the database to default directory, `sphae/workflow/databases`
sphae install

#Install database to specific directory
sphae install --db-dir <directory> 
```

  Install the databases to a directory, `sphae/workflow/databases`

  This workflow requires the 
  - Pfam35.0 database to run viral_verify for contig classification. 
  - CheckV database to test for phage completeness
  - Pharokka databases 
  - Phynteny models
  - Phold databases

This step requires ~17G of storage

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
sphae run --input tests/data/nanopore-subset --preprocess longread --output example --profile slurm -k --threads 16
```

**Command to run only annotation steps**

```bash
#the genomes directory has the already assembled complete genomes
sphae annotate --genome <genomes directory> --output example -k 
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

1. **How do I visualize the phages and gene annotations?**
   To visualize the phages and gene annotations, I recommend using [Clinker](https://github.com/gamcil/clinker). First, gather all the sample genbank files from `sphae.out/RESULTS` and place them in a new directory. Then, execute the clinker command to generate clinker plots, which compare the genes in each genome to each other.
   
   Additionally, for enhanced visualization, consider running [dnaapler](https://github.com/gbouras13/dnaapler) on the genomes in fasta format obtained from 
   `sphae.out/RESULTS`. This step generates reoriented phages that start with terminase genes. Pharokka -> Phold -> Phynteny has to be rerun, and the resulting genbank files can be used for visualization. To perform the annotation steps, run the command 
   `sphae annotate --input <reoriented genomes from dnaapler in fasta format directory>`
   
   Please note that dnaapler may fail if terminase genes are not found, particularly when working with novel phages. The reason these steps haven't been added to sphae. If you encounter any challenges during this process, please feel free to leave an issue, and I'll provide improved documentation to assist you further with the command on how to install and run the command different commands. 

2. **Where are the intermediate files being saved?**
   These files are being saved in `sphae.out/PROCESSING`. If you need more information on the file structure here, or have ideas of better organization then leave an issue and I will make a note to have more documentation. 

3. **Just run annotation on already assembled genomes?**
   
    `sphae annotate --input <input genomes>`
    This command runs only Pharokka, Phold and Phynteny to annotate the assembled genomes. The results are saved to a new directory labeled `sphae.out/annotation`. 

    Note: Currently, Sphae runs Phold in CPU mode, but efforts are underway to support Phold GPU mode for faster processing of this step.

## Issues and Questions

If you come across any issues or errors, report them under [Issues](https://github.com/linsalrob/sphae/issues). 

