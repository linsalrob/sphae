"""

Rules for quality control and quality assurance - Nanopore fastq reads
"""

rule filtlong:
    input:
        i= os.path.join(READDIR, PATTERN)
    output:
        o= os.path.join(QCDIR, "{sample}-filtlong.fastq")
    conda: "../envs/filtlong.yaml"
    log:
        os.path.join(logs, "filtlong_{sample}.log")
    shell:
        """
            export LC_ALL=en_US.UTF-8
            filtlong --min_length 1000 --keep_percent 95 {input.i} > {output.o} 2> {log}
        """

rule rasusa:
    input:
        i= os.path.join(QCDIR, "{sample}-filtlong.fastq")
    output:
        o= os.path.join(QCDIR, "{sample}-rasusa.fastq")
    conda: "../envs/rasusa.yaml"
    log:
        os.path.join(logs, "rasusa_{sample}.log")
    shell:
        """
            rasusa -i {input.i}  --bases 10000000 -o {output.o}
        """