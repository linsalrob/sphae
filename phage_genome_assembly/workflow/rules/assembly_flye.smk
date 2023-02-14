"""

Running flye on nanopore reads 
"""
rule flye:
    input:
        os.path.join(QCDIR, "{sample}-rasusa.fastq")
    threads: 10
    output:
        fasta = os.path.join(ASSEMBLY, "{sample}-flye", "assembly.fasta"),
        gfa = os.path.join(ASSEMBLY, "{sample}-flye", "assembly_graph.gfa"),
        path= os.path.join(ASSEMBLY, "{sample}-flye", "assembly_info.txt")
    params:
        out= os.path.join(ASSEMBLY, "{sample}-flye"),
        model = FLYE_MODEL
    log:
        os.path.join(logs, "flye_{sample}.log")
    conda:
        "../envs/flye.yaml"
    resources:
        mem_mb=64000
    shell:
        """
            flye {params.model} {input} --threads {threads} --out-dir {params.out} 2> {log}
        """