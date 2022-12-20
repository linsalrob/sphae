"""

Running unicycler on nanopore reads separately
"""
rule unicyler_nanopore:
    input:
        s= os.path.join(QCDIR, "{sample}-filtlong.fastq"),
    params:
        out= os.path.join(ASSEMBLY, "{sample}-nanopore-unicycler")
    output:
        fa= os.path.join(ASSEMBLY, "{sample}-nanopore-unicycler", "assembly.fasta"),
        log= os.path.join(ASSEMBLY, "{sample}-nanopore-unicycler", "assembly.gfa")
    log:
        os.path.join(logs, "unicycler_nanopore_{sample}.log")
    threads: 10
    resources:
        mem_mb=64000,
        time=7200,
    conda: "../envs/unicycler.yaml"
    shell:
        """
            unicycler -l {input.s} -o {params.out} -t 10 2> {log}
        """
