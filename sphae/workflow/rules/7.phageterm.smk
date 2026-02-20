"""
Running phageterm on the assembled contigs to identify phage sequences.
"""
rule phageterm_short:
    input:
        r1 = os.path.join(dir_fastp, "{sample}_subsampled_R1.fastq.gz"),
        r2 = os.path.join(dir_fastp, "{sample}_subsampled_R2.fastq.gz"),
        contigs_dir=os.path.join(dir_annot, "{sample}-pr-genomes"),
    output:
        os.path.join(dir_phageterm, "{sample}_pr_phageterm", "{sample}_report.pdf")
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
        cd {params.outdir}
        
        export PYTHONPATH={params.db}:${{PYTHONPATH:-}}
        
        for f in {input.contigs_dir}/*.fasta; do 
            base="$(basename "$f" .fasta)"
                
            phageterm -r "$f" -f {input.r1} -p {input.r2} \
                --report_title {params.outdir}/"$base" \
                -c {threads} 
        done
        touch {output}
        """