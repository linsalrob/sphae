"""
Generating a verial verify search for each contig 
"""
rule viralverify_spades:
    input:
        contigs = os.path.join(dir.spades, "{sample}-pr", "contigs.fasta"),
        db= os.path.join(dir.db, "Pfam35.0", "Pfam-A.hmm.gz")
    output:
        out = os.path.join(dir.spades, "{sample}-pr", "contigs_result_table.csv")
    conda:
        os.path.join(dir.env, "viralverify.yaml")
    params:
        out = os.path.join(dir.spades, "{sample}-pr")
    threads:
        config.resources.smalljob.cpu
    resources:
        mem_mb=config.resources.smalljob.mem,
        time=config.resources.smalljob.time
    log:
        os.path.join(dir.log, "viralverify_spades.{sample}.log")
    benchmark:
        os.path.join(dir.bench, "viralverify_spades.{sample}.txt")
    shell:
        """
        if [[ -s {input.contigs} ]]
        then
            viralverify \
                -f {input.contigs} \
                --hmm {input.db} \
                -o {params.out} \
                -t {threads} \
                &> {log}
        fi
        """


rule viralverify_flye_nano:
    input:
        contigs = os.path.join(dir.flye, "{sample}-sr", "assembly.fasta"),
        db= os.path.join(dir.db, "Pfam35.0", "Pfam-A.hmm.gz")
    output:
        out = os.path.join(dir.flye, "{sample}-sr", "assembly_result_table.csv")
    conda:
        os.path.join(dir.env, "viralverify.yaml")
    params:
        out = os.path.join(dir.flye, "{sample}-sr")
    threads:
        config.resources.smalljob.cpu
    resources:
        mem_mb=config.resources.smalljob.mem,
        time=config.resources.smalljob.time
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
                &> {log}
        fi
        """