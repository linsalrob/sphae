"""
Running phynteny to improve the anntoation 
"""
rule phynteny_run_paired:
    input:
        gbk=os.path.join(dir_annotate, "phold-pr", "{sample}_1_phold", "{sample}_1.gbk"),
    params:
        inputdir=os.path.join(dir_annotate, "{sample}-pr-genomes"),
        idir=os.path.join(dir_annotate, "phold-pr"),
        odir=os.path.join(dir_annotate, "phynteny-pr"),
        model = config['args']['phynteny_db'],
    output:
        pkl=os.path.join(dir_annotate, "phynteny-pr", "{sample}_1_phynteny", "phynteny.gbk")
    conda:
        os.path.join(dir_env, "phynteny.yaml")
    threads:
        config['resources']['smalljob']['cpu']
    resources:
        mem_mb = config['resources']['smalljob']['mem'],
        time = config['resources']['smalljob']['time']
    log:
        os.path.join(dir_log, "phynteny.{sample}.log")
    shell:
        """
        if [[ -s {input.gbk} ]] ; then
            for f in {params.inputdir}/*; do 
                data="$(basename "$f" .fasta)"
                phynteny {params.idir}/"$data"_phold/"$data".gbk -o {params.odir}/"$data"_phynteny \
                    -m {params.model} -f\
                    2> {log}
            done
            touch {output.pkl}
        else
            touch {output.pkl}
        fi
        """

rule phynteny_plotter_paired:
    input:
        gbk=os.path.join(dir_annotate, "phynteny-pr", "{sample}_1_phynteny", "phynteny.gbk"),
        fasta=os.path.join(dir_annotate, "{sample}-pr-genomes" , "{sample}_1.fasta")
    params:
        inputdir=os.path.join(dir_annotate, "{sample}-pr-genomes"),
        idir=os.path.join(dir_annotate, "phynteny-pr"), 
    output:
        plot=os.path.join(dir_annotate, "phynteny-pr", "{sample}_1_phynteny", "plots", "{sample}_1.png")
    resources:
        mem =config['resources']['smalljob']['mem'],
        time = config['resources']['smalljob']['time']
    conda:
        os.path.join(dir_env, "phold.yaml")
    shell:
        """
        if [[ -s {input.gbk} ]] ; then
            for f in {params.inputdir}/*; do 
                data="$(basename "$f" .fasta)"
                genbank_to -g {params.idir}/"$data"_phynteny/phynteny.gbk --gff3 {params.idir}/"$data"_phynteny/phynteny.gff3
                phold plot -i {params.idir}/"$data"_phynteny/phynteny.gbk -f -o {params.idir}/"$data"_phynteny/plots    
            done
            touch {output.plot}
        else
            touch {output.plot}
        fi
        """

rule phynteny_run_nanopore:
    input:
        gbk=os.path.join(dir_annotate, "phold-sr", "{sample}_1_phold", "{sample}_1.gbk"),
    params:
        odir=os.path.join(dir_annotate, "phynteny-sr"),
        idir=os.path.join(dir_annotate, "phold-sr"),
        inputdir=os.path.join(dir_annotate, "{sample}-sr-genomes"),
        model = config['args']['phynteny_db'],
    output:
        pkl=os.path.join(dir_annotate, "phynteny-sr", "{sample}_1_phynteny", "phynteny.gbk")
    conda:
        os.path.join(dir_env, "phynteny.yaml")
    threads:
        config['resources']['smalljob']['cpu']
    resources:
        mem_mb = config['resources']['smalljob']['mem'],
        time = config['resources']['smalljob']['time']
    log:
        os.path.join(dir_log, "phynteny.{sample}.log")
    shell:
        """
        if [[ -s {input.gbk} ]] ; then
            for f in {params.inputdir}/*; do 
                data="$(basename "$f" .fasta)"
                phynteny {params.idir}/"$data"_phold/"$data".gbk -o {params.odir}/"$data"_phynteny \
                    -m {params.model} -f\
                    2> {log}
            done
            touch {output.pkl}
        else
            touch {output.pkl}
        fi
        """

rule phynteny_plotter_longreads:
    input:
        gbk=os.path.join(dir_annotate, "phynteny-sr", "{sample}_1_phynteny", "phynteny.gbk"),
        fasta=os.path.join(dir_annotate, "{sample}-sr-genomes", "{sample}_1.fasta")
    params:
        inputdir=os.path.join(dir_annotate, "{sample}-sr-genomes"),
        idir=os.path.join(dir_annotate, "phynteny-sr"), 
        prefix="phynteny",
    output:
        plot=os.path.join(dir_annotate, "phynteny-sr", "{sample}_1_phynteny", "plots", "{sample}_1.png")
    resources:
        mem =config['resources']['smalljob']['mem'],
        time = config['resources']['smalljob']['time']
    conda:
        os.path.join(dir_env, "phold.yaml")
    shell:
        """
        if [[ -s {input.gbk} ]] ; then
            for f in {params.inputdir}/*; do 
                data="$(basename "$f" .fasta)"
                genbank_to -g {params.idir}/"$data"_phynteny/phynteny.gbk --gff3 {params.idir}/"$data"_phynteny/phynteny.gff3
                phold plot -i {params.idir}/"$data"_phynteny/phynteny.gbk -f -o {params.idir}/"$data"_phynteny/plots                
            done
            touch {output.plot}
        else
            touch {output.plot}
        fi
        """
