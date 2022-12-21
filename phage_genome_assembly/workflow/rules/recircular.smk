"""
Recircularise the phage genomes so they start with terminase large subunit
"""
#rule to map the contig against the terminase subunit
rule terminase_search:
    input:
        contigs = os.path.join(GENOMEDIR, "{sample}.fasta"),
    output:
        tsv = os.path.join(OUTDIR, "recircular", "{sample}-terminase.tsv"),
    params:
        db=os.path.join(DATABASES, 'terminase')
    log:
        os.path.join(logs, "terminase_{sample}.log")
    conda: "../envs/blast.yaml"
    threads: 10
    resources:
        mem_mb=64000
    shell:
        """
            blastx -db {params.db} -query {input} -num_threads {threads} -outfmt 6 -out {output.tsv} 
        """

localrules: rotate_phage
#rotate phage to start with terminase gene 
rule rotate_phage:
    input:
        contigs = os.path.join(GENOMEDIR, "{sample}.fasta"),
        tsv=os.path.join(OUTDIR, "recircular", "{sample}-terminase.tsv"),
    output:
        fa = os.path.join(OUTDIR, "recircular", "{sample}.fasta"),
    log:
        os.path.join(logs, "rotate_{sample}.log")
    conda: "../envs/graph.yaml"
    shell:
        """
            python phage_genome_assembly/workflow/scripts/rotate-phage.py -f {input.contigs} -l {input.tsv} --force > {output.fa} 2> {log}
        """

localrules: rc_phage
#reverse compliment as needed 
rule rc_phage:
    input:
        contigs = os.path.join(OUTDIR, "recircular", "{sample}.fasta"),
    output:
        fa = os.path.join(OUTDIR, "recircular-rc", "{sample}.fasta")
    params:
        indir= os.path.join(OUTDIR, "recircular"),
        outdir= os.path.join(OUTDIR, "recircular-rc")
    log:
        os.path.join(logs, "rc_{sample}.log")
    conda: "../envs/graph.yaml"
    shell:
        """
            python phage_genome_assembly/workflow/scripts/reverse_complement_fasta.py -d {params.indir} -k 8 -o {params.outdir} 2> {log}
        """
