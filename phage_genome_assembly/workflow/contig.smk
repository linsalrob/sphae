"""
phage_genomes
assembly and annotation of phage genomes 

2022, Bhavya Papudeshi

This is an auxiliary Snakefile to install databases or dependencies.
"""


"""CONFIGURATION"""
configfile: os.path.join(workflow.basedir, 'config', 'config.yaml')
configfile: os.path.join(workflow.basedir, 'config', 'databases.yaml')

"""PREFLIGHT CHECKS
Validate your inputs, set up directories, parse your config, etc.
"""
include: "rules/1.preflight-genomes.smk"

"""TARGETS
Declare your targets, either here, or in a separate file.
"""
include: "rules/2.targets-genomes.smk"

"""RULES
Add rules files with the include directive here, or add rules AFTER rule 'all'.
"""

#The rest only if phage genomes have been assembled already
#Recircularise
include: "rules/recircular.smk"

#phage contigs assessment
include: "rules/genome-coverage.smk"
include: "rules/genome-coverage-nano.smk"


"""RUN SNAKEMAKE!"""
rule all:
    input:
        genomes