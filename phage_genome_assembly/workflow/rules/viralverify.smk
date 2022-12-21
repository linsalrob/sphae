"""
Generating a verial verify search for each contig 
"""

rule viralverify_spades:
    input:
        contigs = os.path.join(ASSEMBLY, "{sample}-spades", "contigs.fasta"),
        db= os.path.join(DATABASES, "Pfam35.0", "Pfam-A.hmm.gz")
    output:
        out = os.path.join(ASSEMBLY, "{sample}-viralverify-spades", "contigs_result_table.csv")
    log:
        os.path.join(logs, "viralverify_spades_illumina_{sample}.log")
    conda: "../envs/viralverify.yaml"
    params:
        out = os.path.join(ASSEMBLY, "{sample}-viralverify-spades")
    threads: 10 
    shell:
        """
            if [[ -s {input.contigs} ]]; then
                viralverify -f {input.contigs} --hmm {input.db} -o {params.out} -t {threads} 2> {log}
            fi
        """

rule viralverify_megahit:
    input:
        contigs = os.path.join(ASSEMBLY, "{sample}-megahit", "{sample}.contigs.fa"),
        db= os.path.join(DATABASES, "Pfam35.0", "Pfam-A.hmm.gz")
    output:
        out = os.path.join(ASSEMBLY, "{sample}-viralverify-megahit", "{sample}.contigs_result_table.csv")
    log:
        os.path.join(logs, "viralverify_megahit_illumina_{sample}.log")
    conda: "../envs/viralverify.yaml"
    params:
        out = os.path.join(ASSEMBLY, "{sample}-viralverify-megahit")
    threads: 10 
    shell:
        """
            if [[ -s {input.contigs} ]]; then
                viralverify -f {input.contigs} --hmm {input.db} -o {params.out} -t {threads} 2> {log}
            fi
        """