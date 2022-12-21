"""

Running spades on illumina reads 
"""
rule spades:
    input:
        r1= os.path.join(QCDIR, "{sample}_good_out_R1.fastq"),
        r2= os.path.join(QCDIR, "{sample}_good_out_R2.fastq"),
    params:
        out= os.path.join(ASSEMBLY, "{sample}-spades")
    output:
        fa= os.path.join(ASSEMBLY, "{sample}-spades", "contigs.fasta"),
        gfa= os.path.join(ASSEMBLY, "{sample}-spades", "assembly_graph_with_scaffolds.gfa"),
        path= os.path.join(ASSEMBLY, "{sample}-spades", "contigs.paths")
    log:
        os.path.join(logs, "spades_{sample}.log")
    threads: 10
    resources:
        mem_mb=90000,
        time=7200,
    conda: "../envs/spades.yaml"
    shell:
        """
            spades.py -1 {input.r1} -2 {input.r2} -o {params.out} -t {threads} -m {resources.mem_mb} 2> {log}
        """ 