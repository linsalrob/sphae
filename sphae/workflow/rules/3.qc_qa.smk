"""
Rules for quality control and quality assurance - Illumina paired end reads 
"""
#quality control rules here
rule fastp:
    input:
        r1 = os.path.join(input_dir, PATTERN_R1),
        r2 = os.path.join(input_dir, PATTERN_R2)
    output:
        r1 = os.path.join(dir_fastp,"{sample}_R1.fastq.gz"),
        r2 = os.path.join(dir_fastp,"{sample}_R2.fastq.gz"),
        stats = os.path.join(dir_fastp,"{sample}.stats.json"),
        html = os.path.join(dir_fastp,"{sample}.stats.html")
    conda:
        os.path.join(dir_env, "qc.yaml")
    log:
        os.path.join(dir_log, "fastp.{sample}.log")
    threads: 16
    shell:
        """
        fastp -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} -j {output.stats} -h {output.html} --thread {threads} 2>{log}
        """

rule rasusa:
    input:
        r1 = os.path.join(dir_fastp,"{sample}_R1.fastq.gz"),
        r2 = os.path.join(dir_fastp,"{sample}_R2.fastq.gz"),
    output:
        r1 = os.path.join(dir_fastp,"{sample}_subsampled_R1.fastq.gz"),
        r2 = os.path.join(dir_fastp,"{sample}_subsampled_R2.fastq.gz"),
    params:
        coverage=config['params']['bases'],
        gen_size=config['params']['genomeSize']
    conda:
        os.path.join(dir_env, "qc.yaml")
    log:
        os.path.join(dir_log, "rasusa_paired.{sample}.log")
    shell:
        """
        rasusa reads --coverage {params.coverage} --genome-size {params.gen_size} -o {output.r1} -o {output.r2} {input.r1} {input.r2} 2>{log}
        """


rule filtlong_long:
    """
    runs filtlong to filter quality and length
    """
    input:
        fastq=os.path.join(input_dir, PATTERN_LONG)
    output:
        fastq=os.path.join(dir_nanopore, "{sample}_filt.fastq.gz"),
    conda:
        os.path.join(dir_env, "qc.yaml")
    params:
        qual=config['params']['min_mean_quality'],
        length=config['params']['min_length'],
        target_bases=config['params']['bases']
    log:
        os.path.join(dir_log, "filtlong", "{sample}.log"),
    shell:
        """
        filtlong --target_bases {params.target_bases} --min_mean_q {params.qual} --min_length {params.length} {input.fastq} | pigz > {output.fastq} 2> {log}
        """