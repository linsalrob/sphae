"""
Running CheckV to get the completeness of phage genomes 
"""

rule checkv_megahit:
    input:
        contigs = os.path.join(dir.megahit, "{sample}-pr", "final.contigs.fa")
    output:
        out = os.path.join(dir.megahit, "{sample}-pr", "checkv", "quality_summary.tsv")
    conda:
        os.path.join(dir.env, "checkv.yaml")
    params:
        out = os.path.join(dir.megahit, "{sample}-pr", "checkv"),
        db = os.path.join(dir.db, "checkv-db-v1.5")
    threads:
        config.resources.smalljob.cpu
    resources:
        mem_mb=config.resources.smalljob.mem,
        time=config.resources.smalljob.time
    log:
        os.path.join(dir.log, "checkv_megahit.{sample}.log")
    benchmark:
        os.path.join(dir.bench, "checkv_megahit.{sample}.txt")
    shell:
        """
        export CHECKVDB={params.db}
        if [[ -s {input.contigs} ]] ; then
            checkv end_to_end\
                {input.contigs} \
                {params.out} \
                -t {threads} \
                &> {log}
            touch {output.out}
        else
            touch {output.out}
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
        if [[ -s {input.contigs} ]] ; then
            checkv end_to_end\
                {input.contigs} \
                {params.out} \
                -t {threads} \
                &> {log}
            touch {output.out}
        else
            touch {output.out}
        fi
        """