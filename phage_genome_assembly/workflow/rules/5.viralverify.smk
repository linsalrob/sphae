"""
Generating a verial verify search for each contig 
"""
rule viralverify_megahit:
    input:
        contigs = os.path.join(dir.megahit, "{sample}", "final.contigs.fa"),
        db= os.path.join(dir.db, "Pfam35.0", "Pfam-A.hmm.gz")
    output:
        out = os.path.join(dir.megahit, "{sample}", "final.contigs_result_table.csv")
    conda:
        os.path.join(dir.env, "viralverify.yaml")
    params:
        out = os.path.join(dir.megahit, "{sample}")
    threads:
        config.resources.job.cpu
    resources:
        mem_mb=config.resources.job.mem,
        time=config.resources.job.time
    log:
        os.path.join(dir.log, "viralverify_megahit.{sample}.log")
    benchmark:
        os.path.join(dir.bench, "viralverify_megahit.{sample}.txt")
    shell:
        """
        if [[ -s {input.contigs} ]]
        then
            viralverify \
                -f {input.contigs} \
                --hmm {input.db} \
                -o {params.out} \
                -t {threads} \
                2> {log}
        fi
        """


rule viralverify_flye_nano:
    input:
        contigs = os.path.join(dir.flye, "{sample}", "assembly.fasta"),
        db= os.path.join(dir.db, "Pfam35.0", "Pfam-A.hmm.gz")
    output:
        out = os.path.join(dir.flye, "{sample}", "assembly_result_table.csv")
    conda:
        os.path.join(dir.env, "viralverify.yaml")
    params:
        out = os.path.join(dir.flye, "{sample}")
    threads:
        config.resources.job.cpu
    resources:
        mem_mb=config.resources.job.mem,
        time=config.resources.job.time
    log:
        os.path.join(dir.log, "viralverify_flye_nano.{sample}.log")
    benchmark:
        os.path.join(dir.bench, "viralverify_flye_nano.{sample}.txt")
    shell:
        """
        if [[ -s {input.contigs} ]]
        then
            viralverify \
                -f {input.contigs} \
                --hmm {input.db} \
                -o {params.out} \
                -t {threads} \
                2> {log}
        fi
        """