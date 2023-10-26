"""
Running tax_myphage to genomes
"""
rule taxa_paired:
    """Taxonomic assignment """
    input:
        os.path.join(dir.genome, "{sample}-pr", "{sample}.fasta"),
    output:
        os.path.join(dir.taxa, "{sample}-pr", "Summary_taxonomy.tsv")
    params:
        o=os.path.join(dir.taxa, "{sample}-pr")
    conda:
        os.path.join(dir.env, "taxa.yaml")
    threads:
        config.resources.smalljob.cpu
    resources:
        mem_mb = config.resources.smalljob.mem,
        time = config.resources.smalljob.time
    log:
        os.path.join(dir.log, "taxa.{sample}.log")
    shell:
        """
        taxmyphage -i {input} -t {threads} -o {params.o} --install
        """

rule taxa_flye:
    """Taxonomic assignment """
    input:
        os.path.join(dir.genome, "{sample}-sr", "{sample}.fasta"),
    output:
        os.path.join(dir.taxa, "{sample}-sr", "Summary_taxonomy.tsv")
    params:
        o=os.path.join(dir.taxa, "{sample}-sr")
    conda:
        os.path.join(dir.env, "taxa.yaml")
    threads:
        config.resources.smalljob.cpu
    resources:
        mem_mb = config.resources.smalljob.mem,
        time = config.resources.smalljob.time
    log:
        os.path.join(dir.log, "taxa.{sample}.log")
    shell:
        """
        taxmyphage -i {input} -t {threads} -o {params.o} --install
        """