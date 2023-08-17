"""
Running CheckV to get the completeness of phage genomes 
"""

rule checkv_spades:
    input:
        contigs = os.path.join(dir.spades, "{sample}-pr", "contigs.fasta")
    output:
        out = os.path.join(dir.spades, "{sample}-pr", "checkv", "quality_summary.tsv")
    conda:
        os.path.join(dir.env, "checkv.yaml")
    params:
        out = os.path.join(dir.spades, "{sample}-pr", "checkv"),
        db = os.path.join(dir.db, "checkv-db-v1.5")
    threads:
        config.resources.smalljob.cpu
    resources:
        mem_mb=config.resources.smalljob.mem,
        time=config.resources.smalljob.time
    log:
        os.path.join(dir.log, "checkv_spades.{sample}.log")
    benchmark:
        os.path.join(dir.bench, "checkv_spades.{sample}.txt")
    shell:
        """
        export CHECKVDB={params.db}
        if [[ -s {input.contigs} ]]
        then
            checkv end_to_end\
                {input.contigs} \
                {params.out} \
                -t {threads} \
                &> {log}
        fi
        """


rule checkv_flye_nano:
    input:
        contigs = os.path.join(dir.flye, "{sample}-sr", "assembly.fasta"),
    output:
        out = os.path.join(dir.flye, "{sample}-sr", "checkv", "quality_summary.tsv")
    conda:
        os.path.join(dir.env, "checkv.yaml")
    params:
        out = os.path.join(dir.flye, "{sample}-sr", "checkv"),
        db = os.path.join(dir.db, "checkv-db-v1.5")
    threads:
        config.resources.smalljob.cpu
    resources:
        mem_mb=config.resources.smalljob.mem,
        time=config.resources.smalljob.time
    log:
        os.path.join(dir.log, "checkv_flye_nano.{sample}.log")
    benchmark:
        os.path.join(dir.bench, "checkv_flye_nano.{sample}.txt")
    shell:
        """
        export CHECKVDB={params.db}
        if [[ -s {input.contigs} ]]
        then
            checkv end_to_end\
                {input.contigs} \
                {params.out} \
                -t {threads} \
                &> {log}
        fi
        """