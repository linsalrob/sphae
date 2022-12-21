"""
Generating a table and bam files with read coverage of the phage contig
"""
#rule to build the bam files 
rule genome_coverage_illumina:
    input:
        contigs = os.path.join(OUTDIR, "recircular-rc", "{sample}.fasta"),
        r1= os.path.join(QCDIR, "{sample}_good_out_R1.fastq"),
        r2= os.path.join(QCDIR, "{sample}_good_out_R2.fastq")
    output:
        tsv = os.path.join(OUTDIR, "coverage", "{sample}-illuminaReads.tsv"),
        bam= os.path.join(OUTDIR, "coverage", "{sample}-illuminaReads-bam", "coverm-genome.{sample}_good_out_R1.fastq.bam")
    params:
        bam_dir = os.path.join(OUTDIR, "coverage", "{sample}-illuminaReads-bam")
    log:
        os.path.join(logs, "coverm_ref_illumina_{sample}.log")
    conda: "../envs/coverm.yaml"
    threads: 10
    resources:
        mem_mb=64000
    shell:
        """
            if [[ -s {input.contigs} ]]; then
                export TMPDIR=/scratch/user/nala0006/tmp
                coverm genome -1 {input.r1} -2 {input.r2} --genome-fasta-files {input.contigs} -o {output.tsv} -m coverage_histogram -t {threads} --bam-file-cache-directory {params.bam_dir} 2> {log}
            fi
        """

localrules: index_bam
#indexing the bam files generated 
rule index_bam:
    input:
        bam= os.path.join(OUTDIR, "coverage", "{sample}-illuminaReads-bam", "coverm-genome.{sample}_good_out_R1.fastq.bam")
    output:
        out=os.path.join(OUTDIR, "coverage", "{sample}-illuminaReads-bam", "coverm-genome.{sample}_good_out_R1.fastq.bam.bai")
    conda: "../envs/samtools.yaml"
    shell:
        """
            samtools index {input.bam} {output.out}
        """

#rule to calculate the read coverage per position 
rule read_coverage_illumina:
    input:    
        bam = os.path.join(OUTDIR, "coverage", "{sample}-illuminaReads-bam", "coverm-genome.{sample}_good_out_R1.fastq.bam")
    output:
        tsv=  os.path.join(OUTDIR, "coverage", "{sample}-illuminaReads-bam", "{sample}-Ill-bedtools-genomecov.tsv")
    log:
        os.path.join(logs, "bedtools_Ill_{sample}.log")
    conda: "../envs/bedtools.yaml"
    shell:
        """
            if [[ -s {input.bam} ]]; then
                bedtools genomecov -ibam {input.bam} -d >{output.tsv} 2> {log}
            fi
        """