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
Steps for installing this workflow 

    git clone https://github.com/npbhavya/phage-genomes.git
    cd phage-genomes
    python setup.py install 
    #confirm the workflow is installed by running the below command 
    phage_genome_assembly --help
  
 
## Installing databases
  This workflow requires the 
  - Pfam35.0 database to run viral_verify for contig classification. 
  - large terminase subunit sequences, BLAST database 
  - Pharokka databases 

  This is done automagically when the below command is run. 
  
    phage_genome_assembly install database 
    
## Running the workflow

### Step 1) Assembling phage genoems 

**Steps in this section of workflow**
Pure isolate phages seqeunced on Illumina (paired end) and Nanopore (long read) sequencing technology are processed through the following steps
  - quality control (Illumina: prinseq++, Nanopore: Filtlong)
  - assembly (Illumina: SPAdes and Megahit, Nanopore: Flye and Unicycler)
  - assembly statistics: 
      - read coverage of each contig (CoverM), 
      - contig classification as bacterial/plasmid/viral (viral verify)
      - number of graph components (assembly graph files, python scripts added here)
  
The final ouput is tab separated file providing the summary for each sample assembly, with contig features.

**Commands to run**
Only one command needs to be submitted to run all the above steps: QC, assembly and assembly stats

    #For illumina reads, place the reads both forward and reverse reads to one directory
    phage_genome_assembly run --input test/illumina-subset --output example --profile slurm 

    #For nanopore reads, place the reads, one file per sample in a directory
    phage_genome_assembly run --input test/nanopore-subset --preprocess longread --output example 

    #To run either of the commands on the cluster, add --profile slurm to the command. For instance here is the command for longreads/nanopore reads 
    phage_genome_assembly run --input test/nanopore-subset --preprocess longread --output example --profile slurm 

**Output**

For each sample there should be a tab separated file for each assembler. For instance if test nanopore reads were run through the workflow, then there should be two files within the example/assembly directory
  - reads-assembly-stats_flye.tsv
  - reads-assembly-stats_unicycler.tsv

Each of these files shold contain the 12 columns with the folowing titles and results for each contig assembled. This was a test run, so only one contig was assembled. 

|1    | 2     |    3                                     | 4       |  5       | 6         |  7        |   8      |  9     |  10       |  11 | 12      |
|-----|-------|------------------------------------------|---------|----------|-----------|-----------|----------|--------|-----------|-----|---------|
|Index| Contig| assembly.fasta/reads-filtlong.fastq Mean |Length_x |Circular_x|Connections|Contig name|Prediction|Length_y|Circular_y|Score|Pfam hits|
|0    |contig_1|43.97074                                 |100739   |False     |0          |contig_1   |Virus     |100739  |-         |40.06|Glucosaminidase HNH_3 UDG Asp_protease_2 NUMOD4 GIY-YIG Band_7 Ribonuc_red_lgC DUF1599 Radical_SAM Helicase_C DUF3799 dUTPase Thy1 Toprim_2 NUMOD1 DUF5675 DNA_pol_A_exo1 VWA ThiF AAA_33 NUMOD3 DNA_pol_A |

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
      conda activtate phage-genomes
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

      phage_genome_assembly contig --input test/nanopore-subset --preprocess longread --output ../example --phage-contigs ../example/phage_contigs/ --profile slurm

**Output**

The final output is saved to "example/coverage" directory. This directory includes the following files:
  - tab separated files: includes the number of reads in bases that map to a particular position on the contigs.
  
  | Sample | Genome | coverage | Bases |
  | ---  | ---  |----       | -- |
  | Reads  | Test   | 0        | 0 |
  
  - *.bam drirectory containing two files 
    - bam files
    - *-bedtools-genomecov.tsv
    
     |Genome~Contig name | Genome position | Read coverage|
     | ---- | --- | ---|
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
  
      phage_genome_assembly annotate --phage ../example/phage-final/ --output ../example --profile slurm

**Output**

The output will be saved to "example/pharokka" directory

