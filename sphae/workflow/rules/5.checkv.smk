"""
Running CheckV to get the completeness of phage genomes 
"""

rule checkv_megahit:
    input:
        contigs = os.path.join(dir_megahit, "{sample}-pr", "final.contigs.fa")
    output:
        out = os.path.join(dir_megahit, "{sample}-pr", "checkv", "quality_summary.tsv")
    container:
        "docker://bioedge/checkv:1.0.1"
    params:
        out = os.path.join(dir_megahit, "{sample}-pr", "checkv"),
    threads:
        config['resources']['smalljob']['cpu']
    resources:
        mem_mb=config['resources']['smalljob']['mem'],
        time=config['resources']['smalljob']['time']
    log:
        os.path.join(dir_log, "checkv_megahit.{sample}.log")
    shell:
        """
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
        contigs = os.path.join(dir_flye, "{sample}-sr", "assembly.fasta"),
    output:
        out = os.path.join(dir_flye, "{sample}-sr", "checkv", "quality_summary.tsv")
    container:
        "docker://bioedge/checkv:1.0.1"
    params:
        out = os.path.join(dir_flye, "{sample}-sr", "checkv"),
    threads:
        config['resources']['smalljob']['cpu']
    resources:
        mem_mb=config['resources']['smalljob']['mem'],
        time=config['resources']['smalljob']['time']
    log:
        os.path.join(dir_log, "checkv_flye_nano.{sample}.log")
    shell:
        """
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