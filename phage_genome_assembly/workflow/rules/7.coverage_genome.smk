rule genome_coverage_paired:
    input:
        ref = os.path.join(dir.genome, "{sample}", "{sample}.fasta"),
        reads = os.path.join(dir.temp,"{sample}.paired.tsv")
    output:
        tsv = os.path.join(dir.genome, "{sample}", "megahit_sample_coverage.tsv"),
        bam= os.path.join(dir.genome,"{sample}", "megahit_{sample}.bam"),
        bai= os.path.join(dir.genome,"{sample}", "megahit_{sample}.bam.bai"),
    params:
        out = os.path.join(dir.genome, "{sample}"),
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
        koverage run coverm \
            --reads {input.reads} \
            --ref {input.ref} \
            --output {params.out} \
            --threads {threads} \
            --log {log}
        """


rule genome_coverage_nanopore:
    input:
        ref = os.path.join(dir.genome, "{sample}", "{sample}.fasta"),
        reads = os.path.join(dir.temp, "{sample}.single.tsv")
    output:
        tsv = os.path.join(dir.genome, "{sample}", "flye_sample_coverage.tsv"),
        bam = os.path.join(dir.genome, "{sample}", "flye_{sample}.bam"),
        bai = os.path.join(dir.genome, "{sample}", "flye_{sample}.bam.bai"),
    params:
        out = os.path.join(dir.genome, "{sample}")
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
        koverage run coverm \
            --minimap map-ont \
            --reads {input.reads} \
            --ref {input.ref} \
            --output {params.out} \
            --threads {threads} \
            --log {log}
        """


rule genomecov_paired:
    input:
        bam = os.path.join(dir.genome, "{sample}", "megahit_{sample}.bam"),
        bai = os.path.join(dir.genome, "{sample}", "megahit_{sample}.bam.bai"),
    output:
        tsv=  os.path.join(dir.genome, "{sample}", "megahit_{sample}.gencov.tsv")
    threads:
        config.resources.smalljob.cpu
    resources:
        mem_mb=config.resources.smalljob.mem,
        time=config.resources.smalljob.time
    log:
        os.path.join(dir.log, "genomecov_paired.{sample}.log")
    benchmark:
        os.path.join(dir.bench,"genomecov_paired.{sample}.txt")
    conda:
        os.path.join(dir.env, "bedtools.yaml")
    shell:
        """
        if [[ -s {input.bam} ]]
        then
            bedtools genomecov \
                -ibam {input.bam} \
                -d \
                > {output.tsv} \
                2> {log}
        fi
        """


rule genomecov_nanopore:
    input:
        bam = os.path.join(dir.genome, "{sample}", "flye_{sample}.bam"),
        bai = os.path.join(dir.genome, "{sample}", "flye_{sample}.bam.bai"),
    output:
        tsv=  os.path.join(dir.flye, "{sample}", "flye_{sample}.gencov.tsv")
    threads:
        config.resources.smalljob.cpu
    resources:
        mem_mb=config.resources.smalljob.mem,
        time=config.resources.smalljob.time
    log:
        os.path.join(dir.log, "genomecov_nanopore.{sample}.log")
    benchmark:
        os.path.join(dir.bench,"genomecov_nanopore.{sample}.txt")
    conda:
        os.path.join(dir.env, "bedtools.yaml")
    shell:
        """
        if [[ -s {input.bam} ]]
        then
            bedtools genomecov \
                -ibam {input.bam} \
                -d \
                > {output.tsv} \
                2> {log}
        fi
        """