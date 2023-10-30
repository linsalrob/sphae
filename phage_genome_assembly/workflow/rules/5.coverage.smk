"""
Calculating the read coverage of each contigs
"""
rule write_samples_tsv_paried:
    input:
        r1=lambda wildcards: samples.reads[wildcards.sample]['R1'],
        r2=lambda wildcards: samples.reads[wildcards.sample]['R2']
    output:
        tsv = temp(os.path.join(dir.temp,"{sample}.paired.tsv"))
    localrule: True
    run:
        with open(output.tsv, 'w') as f:
            f.write(f"{wildcards.sample}\t{input.r1}\t{input.r2}\n")


rule write_samples_tsv_single:
    input:
        r1=lambda wildcards: samples.reads[wildcards.sample]['R1'],
    output:
        tsv = temp(os.path.join(dir.temp,"{sample}.single.tsv"))
    localrule: True
    run:
        with open(output.tsv, 'w') as f:
            f.write(f"{wildcards.sample}\t{input.r1}\n")


rule contig_coverage_megahit:
    input:
        ref = os.path.join(dir.megahit, "{sample}-pr", "final.contigs.fa"),
        reads = os.path.join(dir.temp,"{sample}.paired.tsv")
    output:
        tsv = os.path.join(dir.megahit, "{sample}-pr", "results", "sample_coverage.tsv")
    params:
        out = os.path.join(dir.megahit, "{sample}-pr"),
    conda:
        os.path.join(dir.env, "koverage.yaml")
    threads:
        config.resources.smalljob.cpu
    resources:
        mem_mb=config.resources.smalljob.mem,
        time=config.resources.smalljob.time
    log:
        os.path.join(dir.log, "contig_coverage_megahit.{sample}.log")
    benchmark:
        os.path.join(dir.bench, "contig_coverage_megahit.{sample}.txt")
    shell:
        """
        koverage run \
            --reads {input.reads} \
            --ref {input.ref} \
            --output {params.out} \
            --threads {threads} \
            --log {log}
        """


rule contig_coverage_flye_nano:
    input:
        ref = os.path.join(dir.flye,"{sample}-sr", "consensus.fasta"),
        reads = os.path.join(dir.temp, "{sample}.single.tsv")
    output:
        tsv = os.path.join(dir.flye, "{sample}-sr", "results", "sample_coverage.tsv")
    params:
        out = os.path.join(dir.flye, "{sample}-sr")
    conda:
        os.path.join(dir.env, "koverage.yaml")
    threads:
        config.resources.smalljob.cpu
    resources:
        mem_mb=config.resources.smalljob.mem,
        time=config.resources.smalljob.time
    log:
        os.path.join(dir.log, "contig_coverage_flye_nano.{sample}.log")
    benchmark:
        os.path.join(dir.bench, "contig_coverage_flye_nano.{sample}.txt")
    shell:
        """
        koverage run \
            --minimap map-ont \
            --reads {input.reads} \
            --ref {input.ref} \
            --output {params.out} \
            --threads {threads} \
            --log {log}
        """


rule prebuild_koverage:
    output:
        touch(os.path.join(dir.out, "koverage.prebuild"))
    localrule:
        True
    conda:
        os.path.join(dir.env,"koverage.yaml")
    shell:
        """
        koverage run \
            --reads {input.reads} \
            --ref {input.ref} \
            --output {params.out} \
            --threads {threads} \
            --log {log} \
            --conda-create-envs-only
        """