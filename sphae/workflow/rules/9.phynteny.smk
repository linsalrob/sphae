"""
Running phynteny to improve the anntoation 
"""
rule phynteny_run_paired:
    input:
        gbk=os.path.join(dir_pharokka,"{sample}-pr-phold","{sample}.gbk")
    params:
        odir=os.path.join(dir_pharokka, "{sample}-pr-phynteny"),
        model=os.path.join(dir_db, "phynteny_models_zenodo")
    output:
        pkl=os.path.join(dir_pharokka, "{sample}-pr-phynteny", "phynteny.gbk")
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
            phynteny {input.gbk} -o {params.odir} \
                -m {params.model} -f\
                2> {log}
            touch {output.pkl}
        else
            touch {output.pkl}
        fi
        """

rule phynteny_plotter_paired:
    input:
        gbk=os.path.join(dir_pharokka, "{sample}-pr-phynteny", "phynteny.gbk"),
        fasta=os.path.join(dir_genome, "{sample}-pr", "{sample}_genome.fasta")
    params:
        gff3=os.path.join(dir_pharokka, "{sample}-pr-phynteny", "phynteny.gff3"),
        prefix="phynteny",
        output=os.path.join(dir_pharokka, "{sample}-pr-phynteny")
    output:
        plot=os.path.join(dir_pharokka, "{sample}-pr-phynteny", "pharokka_plot.png")
    conda:
        os.path.join(dir_env, "pharokka.yaml")
    shell:
        """
        if [[ -s {input.gbk} ]] ; then
            genbank_to -g {input.gbk} --gff3 {params.gff3}
            pharokka_plotter.py -i {input.fasta} --genbank {input.gbk} --gff {params.gff3} -f -p {params.prefix} -o {params.output}
            touch {output.plot}
        else
            touch {output.plot}
        fi
        """

rule phynteny_run_nanopore:
    input:
        gbk=os.path.join(dir_pharokka,"{sample}-sr-phold","{sample}.gbk")
    params:
        odir=os.path.join(dir_pharokka, "{sample}-sr-phynteny"),
        model=os.path.join(dir_db, "phynteny_models_zenodo")
    output:
        pkl=os.path.join(dir_pharokka, "{sample}-sr-phynteny", "phynteny.gbk")
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
            phynteny {input.gbk} -o {params.odir} \
                -m {params.model} -f\
                2> {log}
            touch {output.pkl}
        else
            touch {output.pkl}
        fi
        """

rule phynteny_plotter_longreads:
    input:
        gbk=os.path.join(dir_pharokka, "{sample}-sr-phynteny", "phynteny.gbk"),
        fasta=os.path.join(dir_genome, "{sample}-sr", "{sample}_genome.fasta")
    params:
        gff3=os.path.join(dir_pharokka, "{sample}-sr-phynteny", "phynteny.gff3"),
        prefix="phynteny",
        output=os.path.join(dir_pharokka, "{sample}-sr-phynteny")
    output:
        plot=os.path.join(dir_pharokka, "{sample}-sr-phynteny", "pharokka_plot.png")
    conda:
        os.path.join(dir_env, "pharokka.yaml")
    shell:
        """
        if [[ -s {input.gbk} ]] ; then
            genbank_to -g {input.gbk} --gff3 {params.gff3}
            pharokka_plotter.py -i {input.fasta} --genbank {input.gbk} --gff {params.gff3} -f -p {params.prefix} -o {params.output}
            touch {output.plot}
        else
            touch {output.plot}
        fi
        """