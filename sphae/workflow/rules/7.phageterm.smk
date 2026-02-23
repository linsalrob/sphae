"""
Running phageterm on the assembled contigs to identify phage sequences.
"""
rule phageterm_short:
    input:
        r1 = os.path.join(dir_fastp, "{sample}_subsampled_R1.fastq.gz"),
        r2 = os.path.join(dir_fastp, "{sample}_subsampled_R2.fastq.gz"),
        contigs_dir=os.path.join(dir_annot, "{sample}-pr-genomes"),
    output:
        os.path.join(dir_phageterm, "{sample}_pr_phageterm", "done.txt")
    params:
        outdir=os.path.join(dir_phageterm, "{sample}_pr_phageterm"),
        db=os.path.join(config['args']['db_dir'], "phagetermvirome-4.3", "phagetermvirome"),
    conda:
        os.path.join(dir_env, "phageterm.yaml")
    threads:
        config['resources']['smalljob']['threads']
    resources:
        mem_mb = config['resources']['smalljob']['mem_mb'],
        runtime = config['resources']['smalljob']['runtime']
    shell:
        """
        mkdir -p {params.outdir}
        export PYTHONPATH={params.db}:${{PYTHONPATH:-}}
        
        for f in {input.contigs_dir}/*; do 
            base="$(basename "$f" .fasta)"

            phageterm -r "$f" -f {input.r1} -p {input.r2} \
                --report_title {params.outdir}/"$base" \
                -c {threads}
            
            mv sphaeoutPROCESSINGph_PhageTerm_report.pdf {params.outdir}/"$base"_report.pdf
            mv sphaeoutPROCESSINGph_sequence.fasta {params.outdir}/"$base"_phageterm.fasta
            mv sphaeoutPROCESSINGph_statistics.csv {params.outdir}/"$base"_phageterm_stats.csv
            pdftotext {params.outdir}/"$base"_report.pdf {params.outdir}/"$base"_report.txt
        done
        touch {output}
        """
