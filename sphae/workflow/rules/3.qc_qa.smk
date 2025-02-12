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
    resources:
        mem =config['resources']['smalljob']['mem'],
        time = config['resources']['smalljob']['time']
    log:
        os.path.join(dir_log, "fastp.{sample}.log")
    threads: 
        config['resources']['smalljob']['cpu']
    shell:
        """
        fastp -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} -j {output.stats} -h {output.html} --thread {threads} 2>{log}
        touch {output.r1}
        touch {output.r2}
        """

rule rasusa:
    input:
        r1 = os.path.join(dir_fastp,"{sample}_R1.fastq.gz"),
        r2 = os.path.join(dir_fastp,"{sample}_R2.fastq.gz"),
    output:
        r1 = os.path.join(dir_fastp,"{sample}_subsampled_R1.fastq.gz"),
        r2 = os.path.join(dir_fastp,"{sample}_subsampled_R2.fastq.gz"),
    params:
        coverage=config['params']['bases']
    conda:
        os.path.join(dir_env, "rasusa.yaml")
    resources:
        mem =config['resources']['smalljob']['mem'],
        time = config['resources']['smalljob']['time']
    threads:
        config['resources']['smalljob']['cpu'],
    log:
        os.path.join(dir_log, "rasusa_paired.{sample}.log")
    shell:
        """
        if [[ -s {input.r1} ]] ; then
            rasusa reads --bases {params.coverage} -o {output.r1} -o {output.r2} {input.r1} {input.r2} 2>{log}
        fi
        touch {output.r1}
        touch {output.r2}
        """

rule run_seqkit_short:
    input:
        r1 = os.path.join(dir_fastp,"{sample}_subsampled_R1.fastq.gz"),
        r2 = os.path.join(dir_fastp,"{sample}_subsampled_R2.fastq.gz"),
    output:
        r=os.path.join(dir_fastp, "{sample}_fastp.txt")
    params:
        r1_temp = os.path.join(dir_fastp, "{sample}_r1.txt"),
        r2_temp = os.path.join(dir_fastp, "{sample}_r2.txt"),
    conda:
        os.path.join(dir_env, "qc.yaml")
    resources:
        mem =config['resources']['smalljob']['mem'],
        time = config['resources']['smalljob']['time']
    threads: 
        config['resources']['smalljob']['cpu'],
    log:
        os.path.join(dir_log, "seqkit", "{sample}_bases_short.log"),
    shell:
        """
        seqkit stats {input.r1} -T > {params.r1_temp}
        seqkit stats {input.r2} -T > {params.r2_temp}

        # Extract numeric values from the second line of the output files
        value1=$(awk 'NR==2 {{print $5}}' {params.r1_temp})
        value2=$(awk 'NR==2 {{print $5}}' {params.r2_temp})

        # Add the values
        sum=$((value1 + value2))

        # Write the sum to output.r
        echo $sum > {output.r}
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
        os.path.join(dir_env, "filtlong.yaml")
    params:
        qual=config['params']['min_mean_quality'],
        length=config['params']['min_length'],
        target_bases=config['params']['bases']
    resources:
        cpu =config['resources']['smalljob']['cpu'],
        mem =config['resources']['smalljob']['mem'],
        time = config['resources']['smalljob']['time']
    log:
        os.path.join(dir_log, "filtlong", "{sample}.log"),
    shell:
        """
        filtlong --target_bases {params.target_bases} --min_mean_q {params.qual} --min_length {params.length} {input.fastq} | pigz > {output.fastq} 2> {log}
        touch {output.fastq}
        """

rule run_seqkit_long:
    input:
        fastq=os.path.join(dir_nanopore, "{sample}_filt.fastq.gz"),
    output:
        r=os.path.join(dir_nanopore, "{sample}_filt.txt"),
    conda:
        os.path.join(dir_env, "qc.yaml")
    resources:
        cpu =config['resources']['smalljob']['cpu'],
        mem =config['resources']['smalljob']['mem'],
        time = config['resources']['smalljob']['time']
    threads:
        config['resources']['smalljob']['cpu'],
    log:
        os.path.join(dir_log, "{sample}_long_seqkit.log"),
    shell:
        """
        seqkit stats {input.fastq} -T -N 50 -N 90 | awk 'NR==2 {{print $5}}' > {output.r}
        """
