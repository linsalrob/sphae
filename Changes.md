# Changes being added

## v1.4.2 
I am updating the default snakemake to v8.11. There have been some changes from v7 to v8, which for now has required me to update 
- how snakemake can submit jobs to the local cluster - https://snakemake.readthedocs.io/en/v8.4.0/executing/cli.html 
- 

## v1.4.1
- added the option to run `sphae annotae` so assembled genomes or reoriented genomes from dnaapler can be run through only the annotation steps
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