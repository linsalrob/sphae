"""

Generating a table with the contig coverages for all the contigs assembled
"""
rule contig_coverage_spades:
    input:
        contigs = os.path.join(ASSEMBLY, "{sample}-spades", "contigs.fasta"),
        r1= os.path.join(QCDIR, "{sample}_good_out_R1.fastq"),
        r2= os.path.join(QCDIR, "{sample}_good_out_R2.fastq")
    output:
        tsv = os.path.join(ASSEMBLY, "{sample}-spades", "{sample}-contigs.tsv")
    log:
        os.path.join(logs, "coverm_spades_illumina_{sample}.log")
    conda: "../envs/coverm.yaml"
    threads: 10
    params:
        tmpdir = TMPDIR
    resources:
        mem_mb=64000
    shell:
        """
            if [[ -s {input.contigs} ]]; then
                export TMPDIR={params.tmpdir}
                coverm contig -1 {input.r1} -2 {input.r2} --reference {input.contigs} -o {output.tsv} -t {threads} 2> {log}
            fi
        """

rule contig_coverage_megahit:
    input:
        contigs = os.path.join(ASSEMBLY, "{sample}-megahit", "{sample}.contigs.fa"),
        r1= os.path.join(QCDIR, "{sample}_good_out_R1.fastq"),
        r2= os.path.join(QCDIR, "{sample}_good_out_R2.fastq")
    output:
        tsv = os.path.join(ASSEMBLY, "{sample}-megahit", "{sample}-contigs.tsv")
    log:
        os.path.join(logs, "coverm_megahit_illumina_{sample}.log")
    conda: "../envs/coverm.yaml"
    threads: 10
    params:
        tmpdir = TMPDIR
    resources:
        mem_mb=64000
    shell:
        """
            if [[ -s {input.contigs} ]]; then
                export TMPDIR={params.tmpdir}
                coverm contig -1 {input.r1} -2 {input.r2} --reference {input.contigs} -o {output.tsv} -t {threads} 2> {log}
            fi
        """