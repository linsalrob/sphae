"""
Running medaka to polish nanopore assembly
"""
rule medaka:
    input:
        fasta = os.path.join(ASSEMBLY, "{sample}-flye", "assembly.fasta"),
        fastq = os.path.join(QCDIR, "{sample}-filtlong.fastq")
    output:
        dir = directory(os.path.join(POLISHING,"{sample}-medaka")),
        fasta = os.path.join(POLISHING,"{sample}-medaka", "consensus.fasta")
    conda:
        "../envs/medaka.yaml"
    params:
        model = MEDAKA_MODEL
    resources:
        mem_mb=32000,
        time=120
    threads:
        16
    shell:
        """
        medaka_consensus -i {input.fastq} -d {input.fasta} -o {output.dir} -m {params.model}  -t {threads}
        """




