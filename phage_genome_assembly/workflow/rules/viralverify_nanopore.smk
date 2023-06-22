"""
Generating a verial verify search for each contig 
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
