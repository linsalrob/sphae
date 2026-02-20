"""
Running phageterm on the assembled contigs to identify phage sequences.
"""
rule phageterm_short:
    input:
        r1 = os.path.join(dir_fastp, "{sample}_subsampled_R1.fastq.gz"),
        r2 = os.path.join(dir_fastp, "{sample}_subsampled_R2.fastq.gz"),
        contigs=os.path.join(dir_annot, "{sample}-pr-genomes", "{sample}_1.fasta"),
    output:
        os.path.join(dir_phageterm, "{sample}_sr_phageterm", "{sample}_report.pdf")
    params:
        inputdir=os.path.join(dir_annot, "{sample}-pr-genomes"),
        outdir=os.path.join(dir_phageterm, "{sample}_sr_phageterm"),
        db=os.path.join(dir_db, "phageterm_db", "phagetermvirome"),
    conda:
        os.path.join(dir_env, "phageterm.yaml")
    threads:
        config['resources']['smalljob']['threads']
    resources:
        mem_mb = config['resources']['smalljob']['mem_mb'],
        runtime = config['resources']['smalljob']['runtime']
    log:
        os.path.join(dir_log, "phageterm_sr.{sample}.log")
    output:
        os.path.join(dir_phageterm, "{sample}", "phageterm_results.tsv")
    shell:
        """
        export PYTHONPATH={params.db}/:$PYTHONPATH
        if [[ -s {input} ]] ; then
            for f in {params.inputdir}/*; do 
                data="$(basename "$f" .fasta)"
                phageterm -r "$f" -f {input.r1} -p {input.r2} \
                    --report_title {params.outdir}/"$f" \
                    -c {threads} > {log} 2>&1
            done
        """