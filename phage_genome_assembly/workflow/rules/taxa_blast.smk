"""
Calculating the blast output against bacteriophage genomes downloaded. 
Note here, the databases were setup manually
Goto: https://www.ncbi.nlm.nih.gov/nuccore
Search for (bacterial virus[Organism]) AND (complete genome[ti]) NOT shotgun[ti] NOT plasmid[ti] NOT bacteria[Organism] 
Dowloaded all the sequences and built a blast database
"""

rule blast:
    input: 
        ins = os.path.join(PHAGEDIR, "{sample}.fasta")
    params:
        #db = os.path.join(DATABASES, "viral_Seq_plusCras")
        db = os.path.join(DATABASES, "viral_Seq")
    output:
        tsv= os.path.join(OUTDIR, "taxa", "{sample}_blastn.tsv")
    threads: 10 
    conda: "../envs/blast.yaml"
    resources:
        mem_mb=64000
    shell:
        """
            blastn -db {params.db} -query {input.ins} -num_threads {threads} -outfmt 6 -out {output.tsv} 
        """

