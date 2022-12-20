"""

Generating a table with the contig coverages for all the contigs assembled - nanopore assemblies 
"""

rule contig_coverage_unicycler_nano:
    input:
        contigs= os.path.join(ASSEMBLY, "{sample}-nanopore-unicycler/assembly.fasta"),
        s= os.path.join(QCDIR, "{sample}-filtlong.fastq")
    output:
        tsv = os.path.join(ASSEMBLY, "{sample}-nanopore-unicycler", "{sample}-contigs.tsv")
    log:
        os.path.join(logs, "coverm_unicycler_nanopore_{sample}.log")
    conda: "../envs/coverm.yaml"
    threads: 10
    resources:
        mem_mb=64000
    shell:
        """
            if [[ -s {input.contigs} ]]; then
                export TMPDIR=/scratch/user/nala0006/tmp
                coverm contig --single {input.s} --reference {input.contigs} -o {output.tsv} -t {threads} 2> {log} 
            fi
        """

rule contig_coverage_flye_nano:
    input:
        contigs = os.path.join(ASSEMBLY, "{sample}-flye/assembly.fasta"),
        s= os.path.join(QCDIR, "{sample}-filtlong.fastq")
    output:
        tsv = os.path.join(ASSEMBLY, "{sample}-flye", "{sample}-contigs.tsv")
    log:
        os.path.join(logs, "coverm_flye_nanopore_{sample}.log")
    conda: "../envs/coverm.yaml"
    threads: 10
    resources:
        mem_mb=64000
    shell:
        """
            if [[ -s {input.contigs} ]]; then
                export TMPDIR=/scratch/user/nala0006/tmp
                coverm contig --single {input.s} --reference {input.contigs} -o {output.tsv} -t {threads} 2> {log}
            fi
        """
