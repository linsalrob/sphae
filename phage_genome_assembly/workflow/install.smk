"""
phage_genomes
assembly and annotation of phage genomes 

2022, Bhavya Papudeshi

This is an auxiliary Snakefile to install databases or dependencies.
"""


import attrmap as ap
import attrmap.utils as au


"""Parse config"""
configfile: os.path.join(workflow.basedir, "..", "config", "config.yaml")
configfile: os.path.join(workflow.basedir, "..", "config", "databases.yaml")
config = ap.AttrMap(config)


"""Preflight"""
include: "rules/1.preflight.smk"
include: "rules/2.targets.smk"


"""RUN SNAKEMAKE"""
rule all:
    input:
        targets.db


"""RULES"""
rule pfam_download:
    params:
        url=os.path.join(config.db.pfam, config.db.pfam_file)
    output:
        os.path.join(dir.db, config.db.pfam_file)
    shell:
        """
            curl -Lo {output} {params.url}
        """


#rule  terminase_download:
#    params:
#        url= os.path.join(config['terminase'])
#    output:
#        o=os.path.join(databaseDir, 'terminase-db2022.zip'),
#    shell:
#        """
#            curl -Lo {output.o} {params.url}
#            unzip {output.o} -d {databaseDir}
#        """


rule  pharokka_download:
    params: 
        pharokka=os.path.join(dir.db, 'pharokka_db')
    output:
        out=os.path.join(dir.db, 'pharokka_db', 'phrogs_db.index')
    conda:
        os.path.join(dir.env, "pharokka.yaml")
    shell:
        """
            install_databases.py -o {params.pharokka}
        """


# rule refseq_mash:
#     params:
#         refseq = os.path.join(dir.db, 'mash_index')
#     output:
#         out=os.path.join(dir.db, 'mash_index', 'refseq.genomes.k21s1000.msh')
#     shell:
#         """
#             wget {params.refseq}
#         """