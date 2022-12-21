"""
Generating a verial verify search for each contig 
"""

rule viralverify_unicycler_nano:
    input:
        db= os.path.join(DATABASES, "Pfam35.0", "Pfam-A.hmm.gz"),
        contig =os.path.join(ASSEMBLY, "{sample}-nanopore-unicycler/assembly.fasta")
    log:
        os.path.join(logs, "viralverify_unicycler_nano_{sample}.log")
    params:
        out = os.path.join(ASSEMBLY, "{sample}-viralverify-unicycler_nano")
    output:
        os.path.join(ASSEMBLY, "{sample}-viralverify-unicycler_nano", "assembly_result_table.csv")
    threads: 10 
    conda: "../envs/viralverify.yaml"
    shell:
        """
            if [[ -s {input.contig} ]]; then
                viralverify -f {input.contig} --hmm {input.db} -o {params.out} -t {threads} 2> {log}
            fi
        """

rule viralverify_flye_nano:
    input:
        contigs = os.path.join(ASSEMBLY, "{sample}-flye/assembly.fasta"),
        db= os.path.join(DATABASES, "Pfam35.0", "Pfam-A.hmm.gz")
    output:
        out = os.path.join(ASSEMBLY, "{sample}-viralverify-flye", "assembly_result_table.csv")
    log:
        os.path.join(logs, "viralverify_flye_nano_{sample}.log")
    conda: "../envs/viralverify.yaml"
    params:
        out = os.path.join(ASSEMBLY, "{sample}-viralverify-flye")
    threads: 10 
    shell:
        """
            if [[ -s {input.contigs} ]]; then
                viralverify -f {input.contigs} --hmm {input.db} -o {params.out} -t {threads} 2> {log}
            fi
        """
