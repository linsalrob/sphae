import yaml
import os
import glob

"""
Parse config
"""
configfile: os.path.join(workflow.basedir, "..", "config", "config.yaml")

"""
DIRECTORIES
"""
dir = {}
dir_env = os.path.join(workflow.basedir, "envs")

# database dir
dir_db = os.path.join(workflow.basedir, 'databases')
print(f"Databases are being saved in, {dir_db} \n")

"""Targets"""
targets = type('', (), {})()
targets.db = []

targets.db.append(os.path.join(dir_db, 'Pfam35.0', 'Pfam-A.hmm.gz'))
targets.db.append(os.path.join(dir_db, 'pharokka_db', 'phrogs_db.index'))
targets.db.append(os.path.join(dir_db, 'checkv-db-v1.5', 'README.txt'))
targets.db.append(os.path.join(dir_db, 'models', 'category_mapping.pkl'))
targets.db.append(os.path.join(dir_db, "phold", "phold_annots.tsv"))
targets.db.append(os.path.join(dir_db, "medaka_models", "medaka.flag"))

"""RUN SNAKEMAKE"""
rule all:
    input:
        targets.db


"""RULES"""
rule pfam_download:
    params:
        os.path.join(config['db']['pfam'], config['db']['pfam_file'])
    output:
        os.path.join(dir_db, config['db']['pfam_file'])
    shell:
        """
            curl -Lo {output} {params}
        """

rule  pharokka_download:
    params: 
        pharokka=os.path.join(dir_db, 'pharokka_db')
    output:
        out=os.path.join(dir_db, 'pharokka_db', 'phrogs_db.index')
    conda:
        os.path.join(dir_env, "pharokka.yaml")
    shell:
        """
            install_databases.py -o {params.pharokka}
        """

rule phynteny_models:
    params:
        models = os.path.join(dir_db)
    output:
        out = os.path.join(dir_db, 'models', 'category_mapping.pkl')
    conda:
        os.path.join(dir_env, "phynteny.yaml")
    shell:
        """
        install_models -o {params.models} -f
        """

rule phold_install:
    params:
        phold_db=os.path.join(dir_db, "phold")
    output:
        os.path.join(dir_db, "phold", "phold_annots.tsv")
    conda:
        os.path.join(dir_env, "phold.yaml")
    shell:
        """
        phold install -d {params.phold_db}
        """

rule checkv_download:
    params:
        os.path.join(dir_db)
    output:
        os.path.join(dir_db, "checkv-db-v1.5", "README.txt")
    conda:
        os.path.join(dir_env, "checkv.yaml")
    shell:
        """
            checkv download_database {params}
        """

rule download_medaka_models:
    params:
        medaka_models=os.path.join(dir_db, "medaka_models")
    output:
        out=os.path.join(dir_db, "medaka_models", "medaka.flag")
    conda:
        os.path.join(dir_env, "medaka.yaml")
    shell:
        """
            medaka tools download_models
            touch {output.out}
        """