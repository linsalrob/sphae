import re
import shutil
import sys
from pathlib import Path

"""
PATTERNS
"""
PATTERN_PROT = "{sample}.faa"

"""
RESOLVER FUNCTION
"""
def resolve_input(wc):
    protein_dir = config['args'].get('proteins')

    protein = os.path.join(protein_dir, f"{wc.sample}.faa") if protein_dir else None

    if protein and Path(protein).exists():
        return protein

    raise ValueError(f"No input found for {wc.sample}")

"""
FUNCTION to get the outputs for pharokka since this is dependent on the input type
"""
rule pharokka_annotate_prot:
    """Annotate genomes with Pharokka for annotate function"""
    input:
        resolve_input,
    params:
        o=os.path.join(dir_annot, "{sample}-pharokka"),
        db = config['args']['pharokka_db'],
        sp="{sample}",
    output:
        faa=os.path.join(dir_annot, "{sample}-pharokka", "{sample}.faa"),
        merged=os.path.join(dir_annot, "{sample}-pharokka", "{sample}_full_merged_output.tsv")
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
            PYTHONWARNINGS="ignore" pharokka_proteins.py \
                -i {input} \
                -o {params.o} \
                -d {params.db} \
                -t {threads} \
                -f -p {params.sp}\
                2> {log}
        fi
        """

rule phold_run_protein:
    input:
        faa=os.path.join(dir_annot, "{sample}-pharokka", "{sample}.faa"),
    params:
        predict=os.path.join(dir_annot, "{sample}-predict"),
        o=os.path.join(dir_annot, "{sample}-phold"),
        prefix="{sample}",
        db = config['args']['phold_db']
    output:
        out=os.path.join(dir_annot, "{sample}-phold","{sample}_aa.fasta"),
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
        if [[ -s {input.faa} ]] ; then
            phold proteins-predict -i {input.faa} -o {params.predict} -p {params.prefix} -t {threads} --cpu -d {params.db} -f 2> {log}
            phold proteins-compare -i {input.faa} --predictions_dir {params.predict} -p {params.prefix} -o {params.o} -t {threads} -d {params.db} -f 2> {log}
        else
            touch {output.out}
        fi
        """