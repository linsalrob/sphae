
rule pharokka:
    """Annotate genomes with Pharokka"""
    input:
        os.path.join(dir.genome, "{sample}.fasta")
    params:
        o=os.path.join(dir.pharokka, "{sample}"),
        db=os.path.join(dir.db, "pharokka_db"),
    output:
        os.path.join(dir.pharokka, "{sample}", "pharokka.gff")
    conda:
        os.path.join(dir.env, "pharokka.yaml")
    threads:
        config.resources.job.cpu
    resources:
        mem_mb = config.resources.job.mem,
        time = config.resources.job.time
    log:
        os.path.join(dir.log, "pharokka.{sample}.log")
    shell:
        """
        pharokka.py \
            -i {input} \
            -o {params.o} \
            -d {params.db} \
            -t {threads} \
            -f \
            2> {log}
        """