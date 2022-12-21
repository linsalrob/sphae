"""

Running flye on nanopore reads 
"""
rule flye:
    input:
        os.path.join(QCDIR, "{sample}-filtlong.fastq")
    threads: 10
    output:
        fasta = os.path.join(ASSEMBLY, "{sample}-flye", "assembly.fasta"),
        gfa = os.path.join(ASSEMBLY, "{sample}-flye", "assembly_graph.gfa"),
        path= os.path.join(ASSEMBLY, "{sample}-flye", "assembly_info.txt")
    params:
        out= os.path.join(ASSEMBLY, "{sample}-flye")
    log:
        os.path.join(logs, "flye_{sample}.log")
    conda:
        "../envs/flye.yaml"
    resources:
        mem_mb=64000
    shell:
        """
            flye --nano-corr {input} --threads {threads} --out-dir {params.out} 2> {log}
        """