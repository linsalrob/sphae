## Miscalleneous scripts

This has a list of misc scripts that can be run on the sphae output 

- num_cds.py : Take a directory of genbank files to count the number of genes that have been assigned to as a gene
  
    `python temp.py -d sphae.output/RESULTS/gbk -o test.csv`

    Output test.csv 
    | Gene  | File1 | File2 | File3 |
    | ------------- | ------------- |
    | DNA primase  | 1 | 1 | 2 |
    | tail protein | 2 | 1 | 1 |
    | holin | 0 | 0 | 0 |


- merging_sphae_output.py: Take a directory of RESULTS from the sphae output, to write the output to a tsv file
  
  For example, 
    `python merging_sphae_output.py KlebPhages_v1.4-out/RESULTS/ test.tsv`

