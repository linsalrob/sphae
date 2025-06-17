# Changes being added

## v1.5.1
- medaka to download models with `sphae install`, now its actually doing this

## v1.5.0
- updated phynteny to phynteny_transformer
- updated medaka to download models with `sphae install`
- updated pharokka yaml file


## v1.4.8
- db_dir variable error fixed 
- moving to using mamba instead of conda as default 
- updating setup.py to pyproject.toml

## v1.4.7 
- checkv database error fixed
- updating phold.yaml and medaka.yaml file 
- updated to set the databases to variables that need to be set by the user
  - this to so the users download only those databases they dont have
  - `sphae install` still works to install all the required databases, but the variables still have to be set to run `sphae run`
    - example of how to set these variables are set in `sphae.sh`  

## v1.4.6
- updating sphae code to make it a smaller container
  - fixing rules to allow --conda-create-envs-only
  - adding a phylogeny module to sphae annotate

## v1.4.5
- catches a specific cases of having mutiple circular phages from assembly
- updating the code to count number of hypothetical proteins to catch other genes that dont have a biological function assigned
- Addressing issue#36- Phrogs annotated toxin not recognised in sphae summary output

## v1.4.4
- adding the option to run pharokka with --pyrodigal-gv to test for alternate coding genes in config file
- sphae plots fix

## v1.4.3
- Summary file update
  - missed adding DTR found or not when only one genome assembled per sample, added this in now
  - If recombinases or transposases are found, the genes are written to the summary. Also if AMR, virulence genes, CRSIPR spacers etc are found. 
    The files are written output. 
- Updating the QC rule to touch the output file so the error is correctly recorded
- the annotate function table after the 3Ps werent being generated, so added that in
- add the number of unknown proteins to the final summary file
- updated phynteny yaml file to include numpy version
- added a python script to the misc folder to merge the RESULTS to one tsv file
- adding the total read length to summary file
- updating the taxa description in the summary file to include the taxa description, lowest taxa and the isolated host from the pharokka inphared result
- added the --no-polish option


## v1.4.2 
I am updating the default snakemake to v8.11. There have been some changes from v7 to v8, which for now has required me to update 
- how snakemake can submit jobs to the local cluster - https://snakemake.readthedocs.io/en/v8.4.0/executing/cli.html 
- The paired end reads have to have the pattern _R1 and _R2 in the sample names, made a note in readme.
- Adding documentation on how to add specific number of base pairs to include in the analysis if phage genome length is known by defining the config file handling multiple phages assembled from one sample
- Checking for DTRs from checkv results to confirm if the contig is an assembled genome.
  
## v1.4.1
- added the option to run `sphae annotate` so assembled genomes or reoriented genomes from dnaapler can be run through only the annotation steps
- Cleaning up code, removed use of Attrmap package
- Updating the summary file to include phage plot from phynteny output
- adding phold results to the summary as needed - making note is phold identified any similar integrase like genes

## v1.4
- Updated the annoation steps to [Pharokka](https://github.com/gbouras13/pharokka), [Phold](https://github.com/gbouras13/phold) and then [Phynteny](https://github.com/susiegriggo/Phynteny)
- Also updated the summary file to confirm the presence/absence of specific genes from the annotation outputs
- Updated Readme

## Upto v1.3.4
Sphae can be run on illumina reads and nanopore reads
- Quality control - [trimnami](https://github.com/beardymcjohnface/Trimnami)
- Assembly - [Megahit](https://github.com/voutcn/megahit) for Illumina reads and [Flye](https://github.com/fenderglass/Flye) + [Medaka](https://github.com/nanoporetech/medaka) for Nanopore reads
  Post assembly the reads are run through [CheckV](https://bitbucket.org/berkeleylab/CheckV/src), [viraverify](https://github.com/ablab/viralVerify) and looking into graph connection to confirm the genome is assembled 
- Annotation - [Pharokka](https://github.com/gbouras13/pharokka) followed by [Phynteny](https://github.com/susiegriggo/Phynteny)
- Generate a summary report file


## Ideas to add
add another module to map the reads to the assembled genome 
  - other fun things to do here
  - calculate the mutations/variants per postion?
