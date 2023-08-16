
rule pharokka_megahit:
    """Annotate genomes with Pharokka"""
    input:
        os.path.join(dir.genome, "{sample}-pr", "{sample}.fasta")
    params:
        o=os.path.join(dir.pharokka, "{sample}-pr"),
        db=os.path.join(dir.db, "pharokka_db"),
    output:
        os.path.join(dir.pharokka, "{sample}-pr", "pharokka.gff")
    conda:
        os.path.join(dir.env, "pharokka.yaml")
    threads:
        config.resources.smalljob.cpu
    resources:
        mem_mb = config.resources.smalljob.mem,
        time = config.resources.smalljob.time
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


rule pharokka_flye:
    """Annotate genomes with Pharokka"""
    input:
        os.path.join(dir.genome, "{sample}-sr", "{sample}.fasta")
    params:
        o=os.path.join(dir.pharokka, "{sample}-sr"),
        db=os.path.join(dir.db, "pharokka_db"),
    output:
        os.path.join(dir.pharokka, "{sample}-sr", "pharokka.gff")
    conda:
        os.path.join(dir.env, "pharokka.yaml")
    threads:
        config.resources.smalljob.cpu
    resources:
        mem_mb = config.resources.smalljob.mem,
        time = config.resources.smalljob.time
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