rule pharokka_annotate:
    """Annotate genomes with Pharokka for annotate function"""
    input:
        os.path.join(input_dir, PATTERN_LONG)
    params:
        o=os.path.join(dir_annot, "{sample}-pharokka"),
        db=os.path.join(dir_db, "pharokka_db"),
        sp="{sample}"
    output:
        gbk=os.path.join(dir_annot, "{sample}-pharokka", "{sample}.gbk"),
        card=os.path.join(dir_annot, "{sample}-pharokka", "top_hits_card.tsv"),
        vfdb=os.path.join(dir_annot, "{sample}-pharokka", "top_hits_vfdb.tsv"),
        spacers=os.path.join(dir_annot, "{sample}-pharokka", "{sample}_minced_spacers.txt"),
        taxa=os.path.join(dir_annot, "{sample}-pharokka", "{sample}_top_hits_mash_inphared.tsv"),
        cdden=os.path.join(dir_annot, "{sample}-pharokka", "{sample}_length_gc_cds_density.tsv"),
        cds=os.path.join(dir_annot, "{sample}-pharokka", "{sample}_cds_functions.tsv")
    conda:
        os.path.join(dir_env, "pharokka.yaml")
    threads:
        config['resources']['smalljob']['cpu']
    resources:
        mem_mb = config['resources']['smalljob']['mem'],
        time = config['resources']['smalljob']['time']
    log:
        os.path.join(dir_log, "pharokka.{sample}.log")
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
            touch {output.gbk}
            touch {output.card}
            touch {output.vfdb}
            touch {output.spacers}
            touch {output.taxa}
            touch {output.cdden}
            touch {output.cds}
        else
            touch {output.gbk}
            touch {output.card}
            touch {output.vfdb}
            touch {output.spacers}
            touch {output.taxa}
            touch {output.cdden}
            touch {output.cds}
        fi
        """

rule phold_run:
    input:
        gbk=os.path.join(dir_annot, "{sample}-pharokka", "{sample}.gbk")
    params:
        predict=os.path.join(dir_annot,"{sample}-predict"),
        o=os.path.join(dir_annot,"{sample}-phold"),
        prefix="{sample}",
        db=os.path.join(dir_db, "phold")
    output:
        gbk=os.path.join(dir_annot,"{sample}-phold","{sample}.gbk"),
        acr=os.path.join(dir_annot,"{sample}-phold","sub_db_tophits", "acr_cds_predictions.tsv"),
        card=os.path.join(dir_annot,"{sample}-phold","sub_db_tophits", "card_cds_predictions.tsv"),
        defense=os.path.join(dir_annot,"{sample}-phold","sub_db_tophits", "defensefinder_cds_predictions.tsv"),
        vfdb=os.path.join(dir_annot,"{sample}-phold","sub_db_tophits", "vfdb_cds_predictions.tsv")
    threads: 
        config['resources']['smalljob']['cpu']
    conda:
        os.path.join(dir_env, "phold.yaml")
    resources:
        mem_mb = config['resources']['smalljob']['mem'],
        time = config['resources']['smalljob']['time']
    log:
        os.path.join(dir_log, "phold.{sample}.log")
    shell:
        """
        if [[ -s {input.gbk} ]] ; then
            phold predict -i {input.gbk} -o {params.predict} -p {params.prefix} -t {threads} --cpu -d {params.db} -f 2> {log}
            phold compare -i {input.gbk} --predictions_dir {params.predict} -p {params.prefix} -o {params.o} -t {threads} -d {params.db} -f 2> {log}
        else
            touch {output.gbk}
        fi
        """

rule phynteny_run:
    input:
        gbk=os.path.join(dir_annot,"{sample}-phold","{sample}.gbk")
    params:
        odir=os.path.join(dir_annot, "{sample}-phynteny"),
        model=os.path.join(dir_db, "phynteny_models_zenodo")
    output:
        pkl=os.path.join(dir_annot, "{sample}-phynteny", "phynteny.gbk")
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

rule phynteny_plotter:
    input:
        gbk=os.path.join(dir_annot, "{sample}-phynteny", "phynteny.gbk"),
        fasta=os.path.join(input_dir, PATTERN_LONG)
    params:
        gff3=os.path.join(dir_annot, "{sample}-phynteny", "phynteny.gff3"),
        prefix="phynteny",
        output=os.path.join(dir_annot, "{sample}-phynteny")
    output:
        plot=os.path.join(dir_annot, "{sample}-phynteny", "pharokka_plot.png")
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

rule summarize_annotations:
    input: 
        pharokka=os.path.join(dir_annot, "{sample}-pharokka", "{sample}.gbk"),
        phold=os.path.join(dir_annot,"{sample}-phold","{sample}.gbk"),
        pkl=os.path.join(dir_annot, "{sample}-phynteny", "phynteny.gbk")
    output:
        pharokka_func=os.path.join(dir_annot, "{sample}-pharokka", "{sample}_pharokka.functions"),
        phold_func=os.path.join(dir_annot,"{sample}-phold","{sample}_phold.functions"),
        pkl_func=os.path.join(dir_annot, "{sample}-phynteny", "phynteny.functions"),
    conda:
        os.path.join(dir_env, "pharokka.yaml")
    shell:
        """
        if [[ -s {input.pharokka} ]] ; then
            genbank_to -g {input.pharokka} -f {output.pharokka_func}
        else
            touch {output.pharokka_func}
        fi

        if [[ -s {input.phold} ]] ; then
            genbank_to -g {input.phold} -f {output.phold_func}
        else
            touch {output.phold_func}
        fi

        if [[ -s {input.pkl} ]] ; then
            genbank_to -g {input.pkl} -f {output.pkl_func}
        else
            touch {output.pkl_func}
        fi
        """

rule annotate_summary:
    input:
        pharokka_func=os.path.join(dir_annot, "{sample}-pharokka", "{sample}_pharokka.functions"),
        phold_func=os.path.join(dir_annot,"{sample}-phold","{sample}_phold.functions"),
        pkl_func=os.path.join(dir_annot, "{sample}-phynteny", "phynteny.functions"),
    output:
        summary_gbk=os.path.join(dir_final, "{sample}", "{sample}_summary.functions")
    params:
        tmp=os.path.join(dir_annot, "{sample}-phynteny", "temp")
    script:
        os.path.join(dir_script, "summary_functions.py")

rule summarize:
    input:
        genome=os.path.join(input_dir, PATTERN_LONG),
        gbk=os.path.join(dir_annot, "{sample}-phynteny", "phynteny.gbk"),
        plot=os.path.join(dir_annot, "{sample}-phynteny", "pharokka_plot.png"),
        ph_taxa =os.path.join(dir_annot, "{sample}-pharokka", "{sample}_top_hits_mash_inphared.tsv"),
        cdden=os.path.join(dir_annot, "{sample}-pharokka", "{sample}_length_gc_cds_density.tsv"),
        cds=os.path.join(dir_annot, "{sample}-pharokka", "{sample}_cds_functions.tsv"),
        amr =os.path.join(dir_annot, "{sample}-pharokka", "top_hits_card.tsv"),
        vfdb=os.path.join(dir_annot, "{sample}-pharokka", "top_hits_vfdb.tsv"),
        spacers=os.path.join(dir_annot, "{sample}-pharokka", "{sample}_minced_spacers.txt"),
        acr=os.path.join(dir_annot,"{sample}-phold","sub_db_tophits", "acr_cds_predictions.tsv"),
        card=os.path.join(dir_annot,"{sample}-phold","sub_db_tophits", "card_cds_predictions.tsv"),
        defense=os.path.join(dir_annot,"{sample}-phold","sub_db_tophits", "defensefinder_cds_predictions.tsv"),
        vfdb_phold=os.path.join(dir_annot,"{sample}-phold","sub_db_tophits", "vfdb_cds_predictions.tsv"),
    output:
        summary=os.path.join(dir_final, "{sample}", "{sample}_summary.txt")
    params:
        genomes= os.path.join(dir_final, "{sample}", "{sample}_genome.fasta"),
        gbks=os.path.join(dir_final, "{sample}", "{sample}.gbk"),
        plots=os.path.join(dir_final, "{sample}", "{sample}_phynteny_plot.png"),
        outdir=os.path.join(dir_final),
        sample="{sample}",
    localrule: True
    script:
        os.path.join(dir_script, 'summary-annot.py')

