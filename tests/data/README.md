# Test dataset for phage_genomes

For Illlumina dataset - 
Downloaded the phix genome, NC_001422.1 Escherichia phage phiX174, complete genome

Then used ART (Artificial Read Generator) using the below command to break the genome to reads

        art_illumina -ss HS25 -l 150 -f 10 -m 200 -s 10 -i genome.fasta -o simulated_reads

The above command, breaks the genome to 150 bp reads, with 200 fragment size, coverage of 10 and standard deviation of 10

This command generated the Illumina test reads, saved within illumina-subset
      


