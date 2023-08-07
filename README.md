# Phage genome toolkit 


**Phage genome assembly and annotation**
This workflow is divided into three sections
1) Assembling the phage isolates using two assemblers and looking at the quality of the assembled contigs \
  Often there are other bacterial, plasmid and prophage contigs along with the phage isolates assembled from the sequences. If there are lots of contigs assembled, then the end of this section requires a manual step. 
  
  Go through the tsv file for each sample, that lists the contig statistics to pick the phage contigs, often these are contigs with highest coverage, look at the genome length - 30 to 50 kbp (or similar to the genome size), and maybe even circular components. 
  Select these contigs and pull them out to new file
  
2) Checking the coverage of the phage genome and variation within the assembled contig
  The selected phage contigs should be saved in fasta format for each sample, and saved in a directory. This section runs coverm again on the contigs to get read coverage across the contig, and recircularises them so they begin with the large terminase gene, and all phages are in the same orientation (no reverse complements).

  Go through the coverage of the phage contigs to confirm there is even coverage. Pick a represenative contig for each sample across assemblies that seems like the highest quality
  
3) Annotation: Running Pharokka 

## Install 
Setting up a new conda environment 

    conda create -n spae 
    conda activate spae

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


## Running the workflow

### Step 1) Assembling phage genomes 

**Steps in this section of workflow**
Isolated phages seqeunced on Illumina (paired end) and Nanopore (long read) sequencing technology are processed through the following steps
  - quality control - trimmnami, https://github.com/beardymcjohnface/Trimnami 
  - assembly (Illumina: SPAdes and Megahit, Nanopore: Flye and Unicycler)
  - assembly statistics: 
      - read coverage of each contig -koverage, https://github.com/beardymcjohnface/Koverage 
      - contig classification as bacterial/plasmid/viral (viral verify)
      - number of graph components (assembly graph files, python scripts added here)
   
The final ouput is tab separated file providing the summary for each sample assembly, with contig features.

**Commands to run**
Only one command needs to be submitted to run all the above steps: QC, assembly and assembly stats

    #For illumina reads, place the reads both forward and reverse reads to one directory
    spae run --input test/illumina-subset --output example

    #For nanopore reads, place the reads, one file per sample in a directory
    spae run --input test/nanopore-subset --preprocess longread --output example 

    #To run either of the commands on the cluster, add --profile slurm to the command. For instance here is the command for longreads/nanopore reads 
    spae run --input test/nanopore-subset --preprocess longread --output example --profile slurm 

**Output**

For each sample there should be a tab separated file for each assembler. For instance if test nanopore reads were run through the workflow, then there should be two files within the example/assembly directory

cd example/assembly
  - reads-assembly-stats_flye.tsv
  - reads-assembly-stats_megahit.tsv

Each of these files shold contain the 20 columns with the folowing titles and results for each contig assembled. This was a test run, so only one contig was assembled. 

Column | Value | Example | Description
--- | --- | --- | ---|
1 | Index | 0 |  |
2 | Sample| reads | Sample name |
2 | Contig | k141_9 | Contig ID from assembly | 
3 | Count | 6070 | Reads mapped to contig |
4 | RPM | 356400.0 | Reads per million |
5 | RPKM | 9610.0 | Reads per kilobase million |
6 | RPK | 163.7 | Reads per kilobase | 
7 | TPM | 41490.0 | Transcripts per million |
8 | Mean | 2714.0 | Estimated mean read depth |
9 | Median | 2798.0 | Estimated median read depth |
10 | Hitrate | 0.9987 | Estimated fraction of contig with depth >0 |
11 | Variance | 1037.0 | Estimated read depth variance |
12 | Length_x | 37091 | contig length |
13 | Circular_x | False | Is the contig circular |
14 | Connections | 31 | Number of unitig connections from assembly graph |
15 | Contig name | k141_9 | Contig ID from assembly | 
16 | Prediction | Virus | Predicted as virus |
17 | Length_y | 37091 | Contig length |
18 | Circular_y | - | Is the contig circular |
19 | Score | 11.96 |  Confidece score of predicting contig as virus, sensitivity threshold | 
20 | Pfam hits | DUF5675 Phage_head_fibr | Genes aligned against Pfam database |
21| contig_id | k141_9 | contig ID from assembly|
22| contig_length |37091 | Contig length |
23| provirus | No| prophage, if there were host boundaries|
24| proviral_length| NA | length of prophage|
25| gene_count| 61| number of genes predicted from this sequence|
26| viral_genes| 6| number of viral genes |
27| host_genes| 0 | number of bacterial sequences|
28| checkv_quality| Medium-quality| classified from the completeness and contamination method|
29| miuvig_quality| genome-fragment | classfied as complete, high-quality or genome-fragment |
30| completeness| 100 | completeness |
31| completeness_method| AAI-based | |
32| contamination| NA | |
33| kmer_freq| | |
34| warnings| | |

### Manual step

**Pick the phage contigs** 

From the output files, pick the phage contigs. These are the contigs that are 
  - Under column "Predictions" or 8, predicted as "Virus"
  - Have the highest coverage: look within column 3
  - Under column 4 or 9, the genome size will be the same as the phage contig length.
  - If circular phage, then under column 5 it would be listed as True.
  
**Separating out only the phage contigs**

Use samtools faidx to grab these contigs 

      #If samtools not installed 
      conda activtate spae
      conda install -c bioconda samtools

      #Run the below command for all the phage contigs per sample
      samtools faidx output/assembly/reads-flye/assembly.fasta 1 >> reads.fasta
      #replace output/assembly/reads-flye/assembly.fasta with a different contigs file 
      #1 is the name of the contig names 
      #reads.fasta, the output fasta file
    
### Step 2) Phage contig quality
The selected phage contigs should be saved in fasta format for each sample, and saved in a directory (example: phage_contigs) 

This section runs coverm again on the contigs to get read coverage across the contig, and recircularises them so they begin with the large terminase gene, and all phages are in the same orientation (no reverse complements).

**Commands to run**

      spae contig --input test/nanopore-subset --preprocess longread --output ../example --phage-contigs ../example/phage_contigs/ --profile slurm

**Output**

The final output is saved to "example/coverage" directory. This directory includes the following files:
  - tab separated files: includes the number of reads in bases that map to a particular position on the contigs.
  
| Sample | Genome | coverage | Bases |
| --- | --- |---| --- |
| Reads | Test | 0        | 0 |
  
- *.bam drirectory containing two files 
  - bam files
  - *-bedtools-genomecov.tsv

|Genome~Contig name | Genome position | Read coverage|
| --- | --- | --- |
 |Test~1 | 1 | 1200 |
     
### Manual step

**Pick a representative contig for each sample with the highest quality** 

From the output files, pick one phage contig per sample. These are the contigs that are 
  - even read coverage across the whole genome 
  - longest contig
  - align the different assemblies and see how many nucleotides are different
  - visualize the read coverge on Tablet
  
  Move the representative assembly from recircular-rc to its own directory (For instance: Phage-contigs-final)
  
### Step 3) Annotation 

**Annotate the phage genomes**

Runs Pharokka to annotate the genomes. 

**Command**

Save the phage genomes to a new directory (in this case, I named the directory phage-final)
  
      spae annotate --phage ../example/phage-final/ --output ../example --profile slurm

**Output**

The output will be saved to "example/pharokka" directory

