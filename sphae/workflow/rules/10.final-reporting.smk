"""
Summarizing the results to one directory
"""
rule summarize_annotations_paired:
    input: 
        pharokka=os.path.join(dir_pharokka, "{sample}-pharokka", "{sample}.gbk"),
        phold=os.path.join(dir_pharokka,"{sample}-pr-phold","{sample}.gbk"),
        pkl=os.path.join(dir_pharokka, "{sample}-pr-phynteny", "phynteny.gbk")
    output:
        pharokka_func=os.path.join(dir_pharokka, "{sample}-pharokka", "{sample}_pr_pharokka.functions"),
        phold_func=os.path.join(dir_pharokka,"{sample}-pr-phold","{sample}_pr_phold.functions"),
        pkl_func=os.path.join(dir_pharokka, "{sample}-pr-phynteny", "phynteny_pr.functions"),
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

rule annotate_summary_paired:
    input:
        pharokka_func=os.path.join(dir_pharokka, "{sample}-pharokka", "{sample}_pr_pharokka.functions"),
        phold_func=os.path.join(dir_pharokka,"{sample}-pr-phold","{sample}_pr_phold.functions"),
        pkl_func=os.path.join(dir_pharokka, "{sample}-pr-phynteny", "phynteny_pr.functions"),
    output:
        summary_gbk=os.path.join(dir_final, "{sample}-pr", "{sample}_summary.functions")
    params:
        tmp=os.path.join(dir_pharokka, "{sample}-pr-phynteny", "temp")
    script:
        os.path.join(dir_script, "summary_functions.py")

rule summarize_paired:
    input:
        assembly=os.path.join(dir_megahit, "{sample}-pr", "log"),
        table= os.path.join(dir_genome, "{sample}-pr", "{sample}-genome-candidates.csv"),
        genome=os.path.join(dir_genome, "{sample}-pr", "{sample}.fasta"),
        gbk=os.path.join(dir_pharokka, "{sample}-pr-phynteny", "phynteny.gbk"),
        plot=os.path.join(dir_pharokka, "{sample}-pr-phynteny", "pharokka_plot.png"),
        ph_taxa =os.path.join(dir_pharokka, "{sample}-pharokka", "{sample}_top_hits_mash_inphared.tsv"),
        cdden=os.path.join(dir_pharokka, "{sample}-pharokka", "{sample}_length_gc_cds_density.tsv"),
        cds=os.path.join(dir_pharokka, "{sample}-pharokka", "{sample}_cds_functions.tsv"),
        amr =os.path.join(dir_pharokka, "{sample}-pharokka", "top_hits_card.tsv"),
        vfdb=os.path.join(dir_pharokka, "{sample}-pharokka", "top_hits_vfdb.tsv"),
        spacers=os.path.join(dir_pharokka, "{sample}-pharokka", "{sample}_minced_spacers.txt"),
        acr=os.path.join(dir_pharokka,"{sample}-pr-phold","sub_db_tophits", "acr_cds_predictions.tsv"),
        card=os.path.join(dir_pharokka,"{sample}-pr-phold","sub_db_tophits", "card_cds_predictions.tsv"),
        defense=os.path.join(dir_pharokka,"{sample}-pr-phold","sub_db_tophits", "defensefinder_cds_predictions.tsv"),
        vfdb_phold=os.path.join(dir_pharokka,"{sample}-pr-phold","sub_db_tophits", "vfdb_cds_predictions.tsv"),
    output:
        summary=os.path.join(dir_final, "{sample}-pr", "{sample}_summary.txt")
    params:
        genomes= os.path.join(dir_final, "{sample}-pr", "{sample}_genome.fasta"),
        gbks=os.path.join(dir_final, "{sample}-pr", "{sample}.gbk"),
        plots=os.path.join(dir_final, "{sample}-pr", "{sample}_phynteny_plot.png"),
        outdir=os.path.join(dir_final),
        sample="{sample}",
    localrule: True
    script:
        os.path.join(dir_script, 'summary.py')

rule summarize_annotations_longreads:
    input: 
        pharokka=os.path.join(dir_pharokka, "{sample}-pharokka", "{sample}.gbk"),
        phold=os.path.join(dir_pharokka,"{sample}-sr-phold","{sample}.gbk"),
        pkl=os.path.join(dir_pharokka, "{sample}-sr-phynteny", "phynteny.gbk")
    output:
        pharokka_func=os.path.join(dir_pharokka, "{sample}-pharokka", "{sample}_sr_pharokka.functions"),
        phold_func=os.path.join(dir_pharokka,"{sample}-sr-phold","{sample}_sr_phold.functions"),
        pkl_func=os.path.join(dir_pharokka, "{sample}-sr-phynteny", "phynteny_sr.functions"),
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

rule annotate_summary_longreads:
    input:
        pharokka_func=os.path.join(dir_pharokka, "{sample}-pharokka", "{sample}_sr_pharokka.functions"),
        phold_func=os.path.join(dir_pharokka,"{sample}-sr-phold","{sample}_sr_phold.functions"),
        pkl_func=os.path.join(dir_pharokka, "{sample}-sr-phynteny", "phynteny_sr.functions"),
    output:
        summary_gbk=os.path.join(dir_final, "{sample}-sr", "{sample}_summary.functions")
    params:
        tmp=os.path.join(dir_pharokka, "{sample}-sr-phynteny", "temp")
    script:
        os.path.join(dir_script, "summary_functions.py")

rule summarize_longread:
    input:
        assembly=os.path.join(dir_flye, "{sample}-sr", "flye.log"),
        table= os.path.join(dir_genome, "{sample}-sr", "{sample}-genome-candidates.csv"),
        genome=os.path.join(dir_genome, "{sample}-sr", "{sample}.fasta"),
        gbk=os.path.join(dir_pharokka, "{sample}-sr-phynteny", "phynteny.gbk"),
        plot=os.path.join(dir_pharokka, "{sample}-sr-phynteny", "pharokka_plot.png"),
        ph_taxa =os.path.join(dir_pharokka, "{sample}-pharokka", "{sample}_top_hits_mash_inphared.tsv"),
        cdden=os.path.join(dir_pharokka, "{sample}-pharokka", "{sample}_length_gc_cds_density.tsv"),
        cds=os.path.join(dir_pharokka, "{sample}-pharokka", "{sample}_cds_functions.tsv"),
        amr =os.path.join(dir_pharokka, "{sample}-pharokka", "top_hits_card.tsv"),
        vfdb=os.path.join(dir_pharokka, "{sample}-pharokka", "top_hits_vfdb.tsv"),
        spacers=os.path.join(dir_pharokka, "{sample}-pharokka", "{sample}_minced_spacers.txt"),
        acr=os.path.join(dir_pharokka,"{sample}-sr-phold","sub_db_tophits", "acr_cds_predictions.tsv"),
        card=os.path.join(dir_pharokka,"{sample}-sr-phold","sub_db_tophits", "card_cds_predictions.tsv"),
        defense=os.path.join(dir_pharokka,"{sample}-sr-phold","sub_db_tophits", "defensefinder_cds_predictions.tsv"),
        vfdb_phold=os.path.join(dir_pharokka,"{sample}-sr-phold","sub_db_tophits", "vfdb_cds_predictions.tsv"),
    output:
        summary=os.path.join(dir_final, "{sample}-sr", "{sample}_summary.txt"),
    params:
        genomes= os.path.join(dir_final, "{sample}-sr", "{sample}_genome.fasta"),
        gbks=os.path.join(dir_final, "{sample}-sr", "{sample}.gbk"),
        sample="{sample}",
        plots=os.path.join(dir_final, "{sample}-sr", "{sample}_phynteny_plot.png"),
        outdir=os.path.join(dir_final),
    localrule: True
    script:
        os.path.join(dir_script, 'summary.py')

