"""

Running megahit on illumina reads
"""
rule megahit:
    input:
        r1= os.path.join(QCDIR, "{sample}_good_out_R1.fastq"),
        r2= os.path.join(QCDIR, "{sample}_good_out_R2.fastq"),
    params:
        out= os.path.join(ASSEMBLY, "{sample}-megahit"),
        n ="{sample}"
    output:
        fa= os.path.join(ASSEMBLY, "{sample}-megahit", "{sample}.contigs.fa")
    log:
        os.path.join(logs, "megahit_{sample}.log")
    threads: 10
    resources:
        mem_mb=90000,
        time=7200,
    conda: "../envs/megahit.yaml"
    shell:
        """
            rmdir {params.out}
            megahit -1 {input.r1} -2 {input.r2} -o {params.out} --out-prefix {params.n} -t {threads} 2> {log}
        """ 

rule fastg:
    input:
        os.path.join(ASSEMBLY, "{sample}-megahit", "{sample}.contigs.fa")
    output:
        os.path.join(ASSEMBLY, "{sample}-megahit", "{sample}.fastg")
    conda: "../envs/megahit.yaml"
    log:
        os.path.join(logs, "megahit_{sample}_fastg.log")
    shell:
        """
            if [[ -s {input} ]]; then
                kmer=$(head -1 {input} | sed 's/>//' | sed 's/_.*//')
                megahit_toolkit contig2fastg $kmer {input} >{output}
            fi 2> {log}
        """