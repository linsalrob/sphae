---
title: ‘Sphae: a bioinformatic workflow to rapidly detect phage therapy candidates among isolated phages’
tags:
  - Python
  - microbiology
  - bacteriophages
  - genomics
  - bioinformatics
authors:
  - name: Bhavya Papudeshi
    orcid: 0000-0001-5359-3100
    corresponding: true
    affiliation: 1 
  - name: George Bouras
     orcid: 0000-0002-5885-4186
    affiliation: 2
  - name: Susanna R. Grigson
     orcid: 0000-0003-4738-3451
     affiliation: 1
  - name: Laura Inglis
    orcid: 0000-0001-7919-8563
    affiliation: 1
  - name: Vijini Mallawaraachchi
     orcid: 0000-0002-2651-8719
     affiliation: 1
 - name: Michael J. Roach
     orcid: 0000-0003-1488-5148
     affiliation: “1, 3”
 - name: Robert A. Edwards
     orcid: 0000-0003-1488-5148
     affiliation: 1 
affiliations:
 - name: Flinders Accelerator for Microbiome Exploration, College of Science and Engineering, Flinders University, Bedford Park, Adelaide, South Australia 5042, Australia   
   index: 1
 - name: Adelaide Medical School, Faculty of Health and Medical Sciences, The University of Adelaide, Adelaide, South Australia 5005, Australia
   index: 2
 - name: Adelaide Centre for Epigenetics and the South Australian Immunogenomics Cancer Institute, Faculty of Health and Medical Sciences, The University  of Adelaide, Adelaide
   index: 3
date: 5 May 2024
bibliography: paper.bib
---
# Summary
The field of phage therapy and antimicrobial research faces significant challenges in characterizing phages, primarily due to the abundance of hypothetical genes and novel genomes. To address this, Sphae offers a streamlined bioinformatic workflow designed for efficient phage characterization and identification of potential therapy candidates. Compatible with Illumina and Oxford Nanopore sequencing platforms, Sphae ensures rigorous quality control and employs advanced assembly and annotation processes. By systematically checking for crucial markers and providing robust error-handling mechanisms, Sphae empowers researchers to make informed decisions about the therapeutic potential of phage isolates. In summary, Sphae represents a significant advancement in phage characterization, facilitating the rapid identification of potential therapeutic interventions amidst the challenges of antimicrobial resistance.
![Sphae workflow.\label{fig:workflow}](Sphae-abstract.png){width=100%}

# Statement of need
Robust tools are urgently needed to aid in the characterisation of bacteriophages in phage therapy and antimicrobial research. With the escalating global challenge of antimicrobial resistance, there is an increasing demand for alternative treatments against bacterial infections. Bacteriophages, owing to their host specificity and ability to circumvent bacterial resistance mechanisms, present a promising avenue for therapy. However, the successful application of phage therapy relies heavily on the precise characterization of phages to identify suitable candidates for treatment. 
A significant challenge in phage characterization lies in the prevalence of phage nomenclature that is essential for effective communication and annotation where several genes are annotated as hypothetical proteins, providing no insight into their functions [@grigson2023knowing]. This lack of detailed genomic information hinders the meticulous selection of appropriate phages for therapeutic use. While sequencing provides comprehensive genomics data, existing bioinformatic tools often lack specificity for the requirement of phage therapy. 
Sphae addresses this need by offering a streamlined bioinformatic workflow tailored for efficient phage characterization and the identification of potential phage therapy candidates. Compatible with both Illumina (paired-end) and Oxford Nanopore (long read) sequencing platforms, Sphae accepts fastq format data as input. Utilizing rigorous quality control methods such as through [Trimnami] (https://github.com/beardymcjohnface/Trimnami), and [Filtlong] (https://github.com/rrwick/Filtlong), Sphae ensures the removal of low-quality reads and adaptor sequences while maintaining data integrity. The resulting reads are subsampled using rasusa [@hall2022rasusa] to include only samples with only 100Mb base pairs from each sample, as 25x to 100x genome coverage is ideal to ensure successful phage assembly [@turner2021phage]. 
The strength of Sphae lies in its assembly and annotation processes. Leveraging tools like MEGAHIT [@li2015megahit] for paired-end read and Flye [@kolmogorov2019assembly] for long read assembly. Although new Nanopore sequencing chemistry is generating higher accuracy, previous chemistry is associated with a relatively higher error rate. To correct for these inaccuracies, Sphae uses [medaka] (https://github.com/nanoporetech/medaka) long read polishing tool. The unique aspect of this workflow is that the assembled contigs are run through [ViralVerify] (https://github.com/ablab/viralVerify) to assign contigs as viral, plasmid or bacteria, CheckV [@nayfach2021checkv] to determine the completeness of the viral contigs and a Python script that calculates contigs connectivity from the assembly graph [@mallawaarachchi2023phables]. These results are combined to determine if the phage genome was assembled, and only the complete phage genome from each sample is moved forward to the annotation step.
Sphae runs Pharokka [@bouras2023pharokka] for comprehensive gene characterisation along with the inbuilt tool to identify the phage’s closely related species in the database. Given that many phage genes are annotated as hypothetical genes, Sphae utilizes [Phold]( https://github.com/gbouras13/phold) that uses a protein language model to generate a 3-D structural alphabet, which is then searched against known structures using FoldSeek [@van2024fast]. The resulting Genbank file with updated annotations is provided as input to Phynteny [@grigson2024phynteny], which uses an LSTM model trained with phage synteny to decipher the functions of annotated genes. Phynteny updates the Genbank further, to reduce the total number of genes annotated as hypothetical genes with an unknown function. From this output, Sphae systematically checks for the presence of crucial markers such as integrase, antimicrobial resistance and virulence genes, essential for determining phage therapy suitability. 
Sphae was built using the Snaketool [@roach2022ten], which streamlines multiple steps into a single bioinformatics workflow. This approach not only enables users to effortlessly download various tools and execute different processes using a single command but also leverages the capabilities of Snakemake to run multiple genomes or workflow steps in parallel. This parallel execution feature proves to be highly efficient, especially when dealing with numerous samples, a trend that is increasingly prevalent in biology research.
```
#Pip install
pip install sphae
#Alternatively, user can use conda or mamba to install 
conda install sphae
#Command to install the databases 
sphae install
#Command to run the program
#For illumina reads
sphae run --input tests/data/illumina-subset --output example -k 
#For nanopore reads
sphae run --input tests/data/nanopore-subset --sequencing longread --output example -k
```
Sphae also uses robust error-handling mechanisms to ensure reliable results, effectively addressing challenges such as incomplete assemblies or fragmented genomes. Through detailed reporting in the final summary file (Table 1), Sphae enables researchers to make informed decisions regarding the therapeutic potential of each phage isolate. 
Table 1: Summary file format
| Attribute                 | Value                         | Attribute explained |
|---------------------------|-------------------------------|--------------|
| Sample                    | nanoreads                     | Sample name |
| Length                    | 100743                        | Length of the phage genome |
| Circular         | False  | Was the genome assembled to be circular, from the assembly graph |
| Graph connections | 0 | In the assembly graph, was the contig connected to other edges. For more information, visualize the assembly graph in GFA format using Bandage| 
| Completeness              | 100.0                 | Phage completeness score from CheckV|
| Contamination             | 0.0                           | Phage contamination score from CheckV|
| Taxa name  | Kehishuvirus sp. 'tikkala'   | Assigned taxa name from Pharokka output | 
| Matching hashes           | 972/1000  | How well did the hashes in the phage genome match to the assigned taxa |
| Number of CDS             | 154                           | Number of genes identified in the genome |
| GC percent                | 0.35                         | 35% GC content
| Coding density | 91.3   | Phages generally have high coding capacity, so if the density is low then something is wrong with the gene calling |
| Integrases    | No  | The presence of integrases could suggest the phage can be integrated into the bacterial host as a prophage |
| Recombinase  | No  | The presence of recombinase could be another indication of the phage being able to maintain as a prophage and also the possibility of homologous and illegitimate recombination|
| Transposase   | No  | The presence of transposase could be an indication of the high rate of gene flow and, the possibility of homologous and illegitimate recombination|
| AMR genes                 | No                            | Antimicrobial resistance genes |
| Virulence factor genes    | No     | Genes that could be playing a role in generating toxins |
| CRISPR spacers | No   | Presence of CRISPR spacers could indicate immunity acquisition to plasmids or other bacteria|

# Availability 
Sphae is distributed on PyPI and as a [Conda](https://conda.io/) package through the Bioconda channel [@Bioconda:2018]. The source code is available on [GitHub](https://github.com/linsalrob/sphae) and features continuous integration tests and continuous deployment using GitHub actions. 

# Acknowledgements
This research/project was undertaken with the assistance of resources and services from Flinders University, Australian Nectar Research Data Commons (ARDC) Nectar Infrastructure, and the National Computational Infrastructure (NCI), supported by the Australian Government. 

# References
Citations to entries in paper.bib 



