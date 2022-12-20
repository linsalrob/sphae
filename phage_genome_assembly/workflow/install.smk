"""
phage_genomes
assembly and annotation of phage genomes 

2022, Bhavya Papudeshi

This is an auxiliary Snakefile to install databases or dependencies.
"""


"""CONFIGURATION"""
configfile: os.path.join(workflow.basedir, 'config', 'config.yaml')
configfile: os.path.join(workflow.basedir, 'config', 'databases.yaml')


include: "rules/1.preflight-database.smk"


"""TARGETS"""
allDatabaseFiles = []

allDatabaseFiles.append(os.path.join(databaseDir, config['pfam_file']))
allDatabaseFiles.append(os.path.join(databaseDir, 'terminase-db2022.zip'))
allDatabaseFiles.append(os.path.join(databaseDir, 'pharokka_db', 'phrogs_db.index'))

"""RUN SNAKEMAKE"""
rule all:
    input:
        allDatabaseFiles


"""RULES"""
rule pfam_download:
    params:
        url=os.path.join(config['pfam'], config['pfam_file'])
    output:
        os.path.join(databaseDir, config['pfam_file'])
    shell:
        """
            curl -Lo {output} {params.url}
        """

rule  terminase_download:
    params:
        url= os.path.join(config['terminase'])
    output:
        o=os.path.join(databaseDir, 'terminase-db2022.zip'),
    shell:
        """
            curl -Lo {output.o} {params.url}
            unzip {output.o} -d {databaseDir}
        """

rule  pharokka_download:
    params: 
        pharokka=os.path.join(databaseDir, 'pharokka_db')
    output:
        out=os.path.join(databaseDir, 'pharokka_db', 'phrogs_db.index')
    conda: "envs/pharokka.yaml"
    shell:
        """
            install_databases.py -o {params.pharokka}
        """

rule refseq_mash:
    params:
        refseq = os.path.join(databaseDir, 'mash_index')
    output:
        out=os.path.join(databaseDir, 'mash_index', 'refseq.genomes.k21s1000.msh')
    shell:
        """
            wget {params.refseq}
        """