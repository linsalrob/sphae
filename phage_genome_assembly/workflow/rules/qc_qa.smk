"""

Rules for quality control and quality assurance - Illumina paired end reads 
"""

rule prinseq:
    input:
        r1 = os.path.join(READDIR, PATTERN_R1),
        r2 = os.path.join(READDIR, PATTERN_R2)
    output:
        r1 = os.path.join(QCDIR, "{sample}_good_out_R1.fastq"),
        r2 = os.path.join(QCDIR, "{sample}_good_out_R2.fastq"),
        s1 = os.path.join(QCDIR, "{sample}_single_out_R1.fastq"),
        s2 = os.path.join(QCDIR, "{sample}_single_out_R2.fastq"),
    conda: "../envs/prinseq.yaml"
    params:
        o = os.path.join(QCDIR, "{sample}")
    log:
        os.path.join(logs, "prinseq_{sample}.log")
    shell:
        """{{
            prinseq++ -min_len 60 -min_qual_mean 25 -ns_max_n 1 -derep 1 \
                    -out_format 0 -trim_tail_left 5 -trim_tail_right 5 \
                    -ns_max_n 5  -trim_qual_type min -trim_qual_left 30 \
                    -trim_qual_right 30 -trim_qual_window 10 \
                    -threads {threads} \
                    -out_name {params.o} \
                    -out_bad /dev/null \
                    -out_bad2 /dev/null \
                    -fastq {input.r1} \
                    -fastq2 {input.r2}; }} 2> {log}
        """
