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