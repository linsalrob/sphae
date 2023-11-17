"""
Rules for quality control and quality assurance - Illumina paired end reads 
"""

rule write_samples_tsv:
    output:
        tsv = temp(os.path.join(dir.temp,"samples.reads.tsv"))
    params:
        sample_dict = samples.reads
    localrule: True
    run:
        from metasnek import fastq_finder
        fastq_finder.write_samples_tsv(params.sample_dict, output.tsv)


rule trimnami:
    input:
        os.path.join(dir.temp,"samples.reads.tsv")
    output:
        targets.qc
    conda:
        os.path.join(dir.env, "trimnami.yaml")
    threads:
        config.resources.bigjob.cpu
    resources:
        mem_mb = config.resources.bigjob.mem,
        time = config.resources.bigjob.time
    params:
        dir = os.path.join(dir.out, "PROCESSING"),
        config = config.args.configfile, 
        trimmer = lambda wildcards: "fastp" if config.args.sequencing == "paired" else "filtlong",
        host = lambda wildcards: "--ref " + config.args.host if config.args.host else "",
        subsample = config.trimnami.qc.subsample,
        profile = lambda wildcards: "--profile " + config.args.profile if config.args.profile else "",
    log:
        os.path.join(dir.log, "trimnami.log")
    benchmark:
        os.path.join(dir.bench,"trimnami_paired.txt")
    shell:
        """
        trimnami run \
            --configfile {params.config} \
            --reads {input} \
            --output {params.dir} \
            --subsample \
            {params.trimmer} \
            {params.host} \
            {params.profile} \
            --log {log}
        """