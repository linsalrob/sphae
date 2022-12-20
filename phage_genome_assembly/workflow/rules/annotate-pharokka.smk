"""

Running pharokka 
"""
rule pharokka:
    input:
        os.path.join(PHAGEDIR, "{sample}.fasta")
    params:
        o=os.path.join(OUTDIR, "pharokka", "{sample}-pharokka"),
        db=os.path.join(DATABASES, "pharokka_db"),
    output:
        os.path.join(OUTDIR, "pharokka", "{sample}-pharokka", "pharokka.gff")
    conda: "../envs/pharokka.yaml"
    log:
        os.path.join(logs, "pharokka_{sample}.log")
    threads: 10
    shell:
        """
            pharokka.py -i {input} -o {params.o} -d {params.db} -t {threads} -f 
        """