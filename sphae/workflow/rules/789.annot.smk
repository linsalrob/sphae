import re
import shutil
import sys
from pathlib import Path
import glob
import os

"""
RESOLVER FUNCTION (USES PREFLIGHT VALUES)
"""

# This will store the sample->file mapping
_sample_to_file = {}

def _init_sample_mapping():
    """Initialize the sample to file mapping from preflight or config."""
    global _sample_to_file
    
    # First, try to import variables from preflight if available
    try:
        # These should be available from 1.preflight-annot.smk since it's included first
        if 'genome_paths' in dir():
            for fp in genome_paths:
                sample_name = os.path.splitext(os.path.basename(fp))[0]
                _sample_to_file[sample_name] = fp
        
        if 'prot_paths' in dir():
            for fp in prot_paths:
                basename = os.path.splitext(os.path.basename(fp))[0]
                sample_name = re.sub(r'[-_]protein$', '', basename)
                if sample_name not in _sample_to_file:
                    _sample_to_file[sample_name] = fp
    except:
        pass
    
    # If still empty, build from config
    if not _sample_to_file:
        genome_dir = config['args'].get('genome')
        protein_dir = config['args'].get('proteins')
        
        if genome_dir:
            genome_paths_local = (
                glob.glob(os.path.join(genome_dir, '*.fasta')) +
                glob.glob(os.path.join(genome_dir, '*.fa')) +
                glob.glob(os.path.join(genome_dir, '*.fna'))
            )
            for fp in genome_paths_local:
                sample_name = os.path.splitext(os.path.basename(fp))[0]
                _sample_to_file[sample_name] = fp
        
        if protein_dir:
            prot_paths_local = glob.glob(os.path.join(protein_dir, '*.faa'))
            for fp in prot_paths_local:
                basename = os.path.splitext(os.path.basename(fp))[0]
                sample_name = re.sub(r'[-_]protein$', '', basename)
                if sample_name not in _sample_to_file:
                    _sample_to_file[sample_name] = fp

def resolve_input_file(wc):
    """Resolve input file for a sample."""
    global _sample_to_file
    
    # Initialize on first call if not already done
    if not _sample_to_file:
        _init_sample_mapping()
    
    if wc.sample in _sample_to_file:
        return _sample_to_file[wc.sample]
    else:
        available = ", ".join(_sample_to_file.keys()) if _sample_to_file else "NONE"
        error_msg = f"\n[ERROR] Cannot find input file for sample: {wc.sample}\n"
        error_msg += f"Available samples: {available}\n"
        error_msg += f"\nSolutions:\n"
        error_msg += f"1. Pass genome directory via command line:\n"
        error_msg += f"   --config genome=/path/to/genomes\n"
        error_msg += f"2. Or update config.yaml with:\n"
        error_msg += f"   args:\n"
        error_msg += f"     genome: /path/to/genomes\n"
        raise FileNotFoundError(error_msg)


def resolve_input_type(wc):
    genome_dir = config['args'].get('genome')
    protein_dir = config['args'].get('proteins')

    if genome_dir:
        if glob.glob(os.path.join(genome_dir, f"{wc.sample}*")):
            return "genome"

    if protein_dir:
        if glob.glob(os.path.join(protein_dir, f"{wc.sample}*")):
            return "protein"

    print(f"[WARNING] No input type found for {wc.sample}")
    return "unknown"

"""
RULES
"""
rule pharokka:
    input:
        infile=resolve_input_file,
        input_type=resolve_input_type
    params:
        o=os.path.join(dir_annot, "{sample}-pharokka"),
        db=config['args']['pharokka_db'],
        sp="{sample}",
        genes=config['params']['gene-predict'],
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
        mem_mb=config['resources']['smalljob']['mem_mb'],
        runtime=config['resources']['smalljob']['runtime']
    log:
        os.path.join(dir_log, "pharokka.{sample}.log")
    run:
        import os
        from pathlib import Path

        infile = input.infile[0] if isinstance(input.infile, list) else input.infile
        input_type = input.input_type

        if not infile or not Path(infile).exists() or Path(infile).stat().st_size == 0:
            shell("""
                touch {output.gbk} {output.card} {output.vfdb} {output.spacers} \
                      {output.taxa} {output.cdden} {output.cds}
            """)
            return

        # decide which mode
        if input_type == "genome":
            cmd = f"""
            PYTHONWARNINGS="ignore" pharokka.py \
                -i {infile} \
                -o {params.o} \
                -d {params.db} \
                -g {params.genes} \
                -t {threads} \
                -f -p {wildcards.sample} \
                2> {log}
            """
        else:
            cmd = f"""
            pharokka_proteins.py \
                -i {infile} \
                -o {params.o} \
                -d {params.db} \
                -t {threads} \
                -f -p {wildcards.sample} \
                2> {log}
            """
        shell(cmd)

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
        gbk=os.path.join(dir_annot, "{sample}-phynteny", "phynteny.gbk")
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
        else
            mkdir -p {output.plot}
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

rule summarize:
    input:
        genome=resolve_input_file,
        input_type=resolve_input_type,
        gbk = os.path.join(dir_annot, "{sample}-phynteny", "phynteny.gbk"),
        plots = os.path.join(dir_annot, "{sample}-phynteny", "plots"),
        ph_taxa = os.path.join(dir_annot, "{sample}-pharokka", "{sample}_top_hits_mash_inphared.tsv"),
        cdden = os.path.join(dir_annot, "{sample}-pharokka", "{sample}_length_gc_cds_density.tsv"),
        cds = os.path.join(dir_annot, "{sample}-pharokka", "{sample}_cds_functions.tsv"),
        amr = os.path.join(dir_annot, "{sample}-pharokka", "top_hits_card.tsv"),
        vfdb = os.path.join(dir_annot, "{sample}-pharokka", "top_hits_vfdb.tsv"),
        spacers = os.path.join(dir_annot, "{sample}-pharokka", "{sample}_minced_spacers.txt"),
        acr = os.path.join(dir_annot, "{sample}-phold","sub_db_tophits", "acr_cds_predictions.tsv"),
        card = os.path.join(dir_annot, "{sample}-phold","sub_db_tophits", "card_cds_predictions.tsv"),
        defense = os.path.join(dir_annot, "{sample}-phold","sub_db_tophits", "defensefinder_cds_predictions.tsv"),
        vfdb_phold = os.path.join(dir_annot, "{sample}-phold","sub_db_tophits", "vfdb_cds_predictions.tsv"),
    output:
        summary = os.path.join(dir_final, "{sample}", "{sample}_summary.txt")
    params:
        genomes = os.path.join(dir_final, "{sample}", "{sample}_input.fasta"),
        gbks = os.path.join(dir_final, "{sample}", "{sample}.gbk"),
        plots = os.path.join(dir_final, "{sample}", "phynteny_plots"),
        sample = "{sample}"
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