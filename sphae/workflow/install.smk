"""
Spae
assembly and annotation of phage genomes 

2022, Bhavya Papudeshi

This is an auxiliary Snakefile to install databases or dependencies.
"""


import attrmap as ap
#import attrmap.utils as au


"""Parse config"""
configfile: os.path.join(workflow.basedir, "..", "config", "config.yaml")
config = ap.AttrMap(config)



"""Preflight"""
include: "rules/1.preflight-database.smk"

"""Targets"""
targets = ap.AttrMap()

targets.db = []

targets.db.append(os.path.join(databaseDir, 'Pfam35.0', 'Pfam-A.hmm.gz'))
targets.db.append(os.path.join(databaseDir, 'pharokka_db', 'phrogs_db.index'))
targets.db.append(os.path.join(databaseDir, 'checkv-db-v1.5', 'README.txt'))
targets.db.append(os.path.join(databaseDir, 'phynteny_models_zenodo', 'grid_search_model.m_400.b_256.lr_0.0001.dr_0.1.l_2.a_tanh.o_rmsprop.rep_0.best_val_loss.h5'))


"""RUN SNAKEMAKE"""
rule all:
    input:
        targets.db


"""RULES"""
rule pfam_download:
    params:
        os.path.join(config.db.pfam, config.db.pfam_file)
    output:
        os.path.join(databaseDir, config.db.pfam_file)
    shell:
        """
            curl -Lo {output} {params}
        """

rule  pharokka_download:
    params: 
        pharokka=os.path.join(databaseDir, 'pharokka_db')
    output:
        out=os.path.join(databaseDir, 'pharokka_db', 'phrogs_db.index')
    conda:
        os.path.join(dir.env, "pharokka.yaml")
    shell:
        """
            install_databases.py -o {params.pharokka}
        """

rule checkv_database:
    params:
        checkv_db=os.path.join(databaseDir)
    output:
        os.path.join(databaseDir, "checkv-db-v1.5", "README.txt")
    conda:
        os.path.join(dir.env, "checkv.yaml")
    shell:
        """
            checkv download_database {params.checkv_db}
        """

rule phynteny_models:
    params:
        url = os.path.join(config.db.models),
        download = os.path.join(databaseDir, "phynteny_models_v0.1.11.tar.gz"),
        models = os.path.join(databaseDir)
    output:
        out = os.path.join(databaseDir, 'phynteny_models_zenodo', 'grid_search_model.m_400.b_256.lr_0.0001.dr_0.1.l_2.a_tanh.o_rmsprop.rep_0.best_val_loss.h5')
    conda:
        os.path.join(dir.env, "phynteny.yaml")
    shell:
        """
            wget -O {params.download} {params.url}
            tar -xvzf {params.download} -C {params.models}
        """
