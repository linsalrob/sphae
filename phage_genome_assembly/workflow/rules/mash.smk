"""
Calculating the mash distances of the phage genome asssembled against 
the refeseq genomes
"""

rule mash_index:
    input:
        os.path.join(PHAGEDIR, "{sample}.fasta")
    log:
        os.path.join(logs, "mash_index_{sample}.log")
    output:
        os.path.join(OUTDIR, "taxa", "{sample}.msh")
    threads: 10 
    conda: "../envs/mash.yaml"
    shell:
        """
            mash sketch -m 2 {input} -o {output}
        """

rule mash_distance:
    input:
        os.path.join(OUTDIR, "taxa", "{sample}.msh")
    params:
        db=os.path.join(DATABASES, "mash_index", "refseq.genomes.k21s1000.msh")
    output:
        os.path.join(OUTDIR, "taxa", "{sample}_distances.tab")
    conda: "../envs/mash.yaml"
    log:
        os.path.join(logs, "mash_distance_{sample}.log")
    shell:
        """
            mash dist {params.db} {input} > {output}
        """