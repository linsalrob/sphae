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
try:
    if config['args']['db_dir'] is None:
        dir_db = os.path.join(workflow.basedir, 'databases')
    else:
        dir_db = config['args']['db_dir']
except KeyError:
    dir_db = os.path.join(workflow.basedir, 'databases')
print(f"Databases are being saved in, {dir_db} \n")

"""Targets"""
targets = type('', (), {})()
targets.db = []

targets.db.append(os.path.join(dir_db, 'Pfam35.0', 'Pfam-A.hmm.gz'))
targets.db.append(os.path.join(dir_db, 'pharokka_db', 'phrogs_db.index'))
targets.db.append(os.path.join(dir_db, 'checkv-db-v1.5', 'README.txt'))
targets.db.append(os.path.join(dir_db, 'phynteny_models_zenodo', 'grid_search_model.m_400.b_256.lr_0.0001.dr_0.1.l_2.a_tanh.o_rmsprop.rep_0.best_val_loss.h5'))
targets.db.append(os.path.join(dir_db, "phold", "phold_annots.tsv"))

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

rule checkv_database:
    params:
        checkv_db=os.path.join(dir_db)
    output:
        os.path.join(dir_db, "checkv-db-v1.5", "README.txt")
    conda:
        os.path.join(dir_env, "checkv.yaml")
    shell:
        """
            checkv download_database {params.checkv_db}
        """

rule phynteny_models:
    params:
        url = os.path.join(config['db']['models']),
        download = os.path.join(dir_db, "phynteny_models_v0.1.11.tar.gz"),
        models = os.path.join(dir_db)
    output:
        out = os.path.join(dir_db, 'phynteny_models_zenodo', 'grid_search_model.m_400.b_256.lr_0.0001.dr_0.1.l_2.a_tanh.o_rmsprop.rep_0.best_val_loss.h5')
    conda:
        os.path.join(dir_env, "phynteny.yaml")
    shell:
        """
            wget -O {params.download} {params.url}
            tar -xvzf {params.download} -C {params.models}
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

