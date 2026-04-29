import re
import shutil
import sys
from pathlib import Path

"""
PATTERNS
"""
PATTERN_LONG = "{sample}.fasta"
PATTERN_PROT = "{sample}.faa"

"""
RESOLVER FUNCTION
"""
def resolve_input(wc):
    genome_dir = config['args'].get('genome')
    protein_dir = config['args'].get('proteins')

    genome = os.path.join(genome_dir, f"{wc.sample}.fasta") if genome_dir else None
    protein = os.path.join(protein_dir, f"{wc.sample}.faa") if protein_dir else None

    if genome and Path(genome).exists():
        return genome
    if protein and Path(protein).exists():
        return protein

    raise ValueError(f"No input found for {wc.sample}")

if config['args'].get('genome'):
    input_type = "genome"
elif config['args'].get('proteins'):
    input_type = "proteins"
else:
    raise ValueError("Please provide either a genome directory or a protein directory, but not both.")

rule pharokka_annotate_genome:
    """Annotate genomes with Pharokka for annotate function"""
    input:
        resolve_input
    params:
        itype=input_type,
        o=os.path.join(dir_annot, "{sample}-pharokka"),
        db = config['args']['pharokka_db'],
        sp="{sample}",
        genes= config['params']['gene-predict'],
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
        config['resources']['smalljob']['threads']
    resources:
        mem_mb = config['resources']['smalljob']['mem_mb'],
        runtime = config['resources']['smalljob']['runtime']
    log:
        os.path.join(dir_log, "pharokka.{sample}.log")
    shell:
        """
        if [[ -s {input} ]] ; then
            if [[ "{params.itype}" == "genome" ]]; then
                PYTHONWARNINGS="ignore" pharokka.py \
                    -i {input} \
                    -o {params.o} \
                    -d {params.db} \
                    -g {params.genes} \
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
                PYTHONWARNINGS="ignore" pharokka_proteins.py \
                    -i {input} \
                    -o {params.o} \
                    -d {params.db} \
                    -t {threads} \
                    -f -p {params.sp} \
                    2> {log}
            fi
        fi
        """

rule phold_run:
    input:
        gbk=os.path.join(dir_annot, "{sample}-pharokka", "{sample}.gbk")
    params:
        predict=os.path.join(dir_annot, "{sample}-predict"),
        o=os.path.join(dir_annot, "{sample}-phold"),
        prefix="{sample}",
        db = config['args']['phold_db']
    output:
        gbk=os.path.join(dir_annot, "{sample}-phold","{sample}.gbk"),
        acr=os.path.join(dir_annot, "{sample}-phold","sub_db_tophits", "acr_cds_predictions.tsv"),
        card=os.path.join(dir_annot, "{sample}-phold","sub_db_tophits", "card_cds_predictions.tsv"),
        defense=os.path.join(dir_annot, "{sample}-phold","sub_db_tophits", "defensefinder_cds_predictions.tsv"),
        vfdb=os.path.join(dir_annot, "{sample}-phold","sub_db_tophits", "vfdb_cds_predictions.tsv")
    threads: 
        config['resources']['smalljob']['threads']
    conda:
        os.path.join(dir_env, "phold.yaml")
    resources:
        mem_mb = config['resources']['smalljob']['mem_mb'],
        runtime = config['resources']['smalljob']['runtime']
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
        gbk=os.path.join(dir_annot, "{sample}-phold","{sample}.gbk")
    params:
        odir=os.path.join(dir_annot, "{sample}-phynteny"),
        model = config['args']['phynteny_db'],
    output:
        pkl=os.path.join(dir_annot, "{sample}-phynteny", "phynteny.gbk")
    conda:
        os.path.join(dir_env, "phynteny.yaml")
    threads:
        config['resources']['smalljob']['threads']
    resources:
        mem_mb = config['resources']['smalljob']['mem_mb'],
        runtime = config['resources']['smalljob']['runtime']
    log:
        os.path.join(dir_log, "phynteny.{sample}.log")
    shell:
        """
        if [[ -s {input.gbk} ]] ; then
            phynteny_transformer  {input.gbk} -o {params.odir} \
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
        fasta=resolve_input
    params:
        gff3=os.path.join(dir_annot, "{sample}-phynteny", "phynteny.gff3"),
        prefix="{sample}",
        output=os.path.join(dir_annot, "{sample}-phynteny", "plots")
    output:
        plot=directory(os.path.join(dir_annot, "{sample}-phynteny", "plots"))
    resources:
        mem_mb =config['resources']['smalljob']['mem_mb'],
        runtime = config['resources']['smalljob']['runtime']
    conda:
        os.path.join(dir_env, "phold.yaml")
    shell:
        """
        if [[ -s {input.gbk} ]] ; then
            genbank_to -g {input.gbk} --gff3 {params.gff3}
            phold plot -i {input.gbk} -f -p {params.prefix} -o {params.output}
        fi
        """

rule summarize_annotations:
    input: 
        pharokka=os.path.join(dir_annot, "{sample}-pharokka", "{sample}.gbk"),
        phold=os.path.join(dir_annot, "{sample}-phold","{sample}.gbk"),
        pkl=os.path.join(dir_annot, "{sample}-phynteny", "phynteny.gbk")
    output:
        pharokka_func=os.path.join(dir_annot, "{sample}-pharokka", "{sample}_pharokka.functions"),
        phold_func=os.path.join(dir_annot, "{sample}-phold","{sample}_phold.functions"),
        pkl_func=os.path.join(dir_annot, "{sample}-phynteny", "phynteny.functions"),
    resources:
        mem_mb =config['resources']['smalljob']['mem_mb'],
        runtime = config['resources']['smalljob']['runtime']
    conda:
        os.path.join(dir_env, "phold.yaml")
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
    localrule: True
    script:
        os.path.join(dir_script, "summary_annot_functions.py")

rule summarize:
    input:
        genome=resolve_input,
        gbk=os.path.join(dir_annot, "{sample}-phynteny", "phynteny.gbk"),
        plot=directory(os.path.join(dir_annot, "{sample}-phynteny", "plots")),
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
        plots=directory(os.path.join(dir_final, "{sample}", "{sample}_phynteny")),
        outdir=os.path.join(dir_final),
        sample="{sample}",
    localrule: True
    script:
        os.path.join(dir_script, 'summary-annot.py')


rule accessory_files:
    input:
        summary=os.path.join(dir_final, "{sample}", "{sample}_summary.txt"),
        acr=os.path.join(dir_annot,"{sample}-phold","sub_db_tophits", "acr_cds_predictions.tsv"),
        card=os.path.join(dir_annot,"{sample}-phold","sub_db_tophits", "card_cds_predictions.tsv"),
        defense=os.path.join(dir_annot,"{sample}-phold","sub_db_tophits", "defensefinder_cds_predictions.tsv"),
        vfdb_phold=os.path.join(dir_annot,"{sample}-phold","sub_db_tophits", "vfdb_cds_predictions.tsv"),
    output:
        amr=os.path.join(dir_final, "{sample}", "{sample}_phold_amr.tsv"),
        vfdb=os.path.join(dir_final, "{sample}", "{sample}_phold_vfdb.tsv"),
        acr=os.path.join(dir_final, "{sample}", "{sample}_phold_acr.tsv"),
        defense=os.path.join(dir_final, "{sample}", "{sample}_phold_defense.tsv"),
    params:
        sample="{sample}"
    localrule: True
    shell:
        """
        cp {input.acr} {output.acr}
        cp {input.card} {output.amr}
        cp {input.defense} {output.defense}
        cp {input.vfdb_phold} {output.vfdb}
        """