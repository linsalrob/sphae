rule pharokka_annotate:
    """Annotate genomes with Pharokka for annotate function"""
    input:
        os.path.join(samples.reads, "{sample}.fasta")
    params:
        o=os.path.join(dir.pharokka, "{sample}"),
        db=os.path.join(dir.db, "pharokka_db"),
        sp="{sample}"
    output:
        gbk=os.path.join(dir.pharokka, "{sample}", "{sample}.gbk"),
        plot=os.path.join(dir.pharokka, "{sample}", "{sample}_pharokka_plot.png"),
        card=os.path.join(dir.pharokka, "{sample}", "top_hits_card.tsv"),
        vfdb=os.path.join(dir.pharokka, "{sample}", "top_hits_vfdb.tsv"),
        spacers=os.path.join(dir.pharokka, "{sample}", "{sample}_minced_spacers.txt"),
        taxa=os.path.join(dir.pharokka, "{sample}", "{sample}_top_hits_mash_inphared.tsv"),
        cdden=os.path.join(dir.pharokka, "{sample}", "{sample}_length_gc_cds_density.tsv"),
        cds=os.path.join(dir.pharokka, "{sample}", "{sample}_cds_functions.tsv")
    conda:
        os.path.join(dir.env, "pharokka.yaml")
    threads:
        config.resources.smalljob.cpu
    resources:
        mem_mb = config.resources.smalljob.mem,
        time = config.resources.smalljob.time
    log:
        os.path.join(dir.log, "pharokka_annotate.{sample}.log")
    benchmark:
        os.path.join(dir.bench,"pharokka_annotate_{sample}.txt")
    shell:
        """
        if [[ -s {input} ]] ; then
            pharokka.py \
                -i {input} \
                -o {params.o} \
                -d {params.db} \
                -t {threads} \
                -f -p {params.sp}\
                2> {log}

            pharokka_plotter.py -i {input} -n {params.sp}_pharokka_plot -o {params.o} -p {params.sp} -f
            touch {output.gbk}
            touch {output.plot}
            touch {output.card}
            touch {output.vfdb}
            touch {output.spacers}
            touch {output.taxa}
            touch {output.cdden}
            touch {output.cds}
        else
            touch {output.gbk}
            touch {output.plot}
            touch {output.card}
            touch {output.vfdb}
            touch {output.spacers}
            touch {output.taxa}
            touch {output.cdden}
            touch {output.cds}
        fi
        """

rule phold_annotate:
    input:
        gbk=os.path.join(dir.pharokka, "{sample}", "{sample}.gbk")
    params:
        predict=os.path.join(dir.pharokka,"{sample}-predict"),
        o=os.path.join(dir.pharokka,"{sample}-phold"),
        prefix="{sample}",
        db=os.path.join(dir.db, "phold")
    output:
        gbk=os.path.join(dir.pharokka,"{sample}-phold","{sample}.gbk"),
    threads: 
        config.resources.smalljob.cpu
    conda:
        os.path.join(dir.env, "phold.yaml")
    resources:
        mem_mb = config.resources.smalljob.mem,
        time = config.resources.smalljob.time
    log:
        os.path.join(dir.log, "phold.{sample}.log")
    benchmark:
        os.path.join(dir.bench,"phold_{sample}.txt")
    shell:
        """
        if [[ -s {input.gbk} ]] ; then
            phold predict -i {input.gbk} -o {params.predict} -p {params.prefix} -t {threads} --cpu -d {params.db} -f
            phold compare -i {input.gbk} --predictions_dir {params.predict} -p {params.prefix} -o {params.o} -t {threads} -d {params.db} -f
        else
            touch {output.gbk}
        fi
        """
    
rule phynteny_annotate:
    input:
        gbk=os.path.join(dir.pharokka,"{sample}-phold","{sample}.gbk")
    params:
        odir=os.path.join(dir.pharokka, "{sample}", "phynteny"),
        model=os.path.join(dir.db, "phynteny_models_zenodo")
    output:
        pkl=os.path.join(dir.pharokka, "{sample}", "phynteny", "phynteny.gbk")
    conda:
        os.path.join(dir.env, "phynteny.yaml")
    threads:
        config.resources.smalljob.cpu
    resources:
        mem_mb = config.resources.smalljob.mem,
        time = config.resources.smalljob.time
    log:
        os.path.join(dir.log, "phynteny.{sample}.log")
    benchmark:
        os.path.join(dir.bench,"phynteny_nanopore_{sample}.txt")
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