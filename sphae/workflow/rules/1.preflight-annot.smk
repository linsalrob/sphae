import os
import glob
import yaml
import re
import shutil
import sys
from pathlib import Path

"""
CONFIG FILE
"""
configfile: os.path.join(workflow.basedir, "..", "config", "config.yaml")

"""
DIRECTORIES
"""
dir = {}

"""
INPUT DIRECTORIES
"""
genome_dir = config['args'].get('genome')
protein_dir = config['args'].get('proteins')  # add this in config

"""
READ GENOME FILES
"""
genome_paths = []
samples_genome = []

if genome_dir:
    genome_paths = (
        glob.glob(os.path.join(genome_dir, '*.fasta')) +
        glob.glob(os.path.join(genome_dir, '*.fa')) +
        glob.glob(os.path.join(genome_dir, '*.fna'))
    )

    samples_genome = [
        os.path.splitext(os.path.basename(fp))[0]
        for fp in genome_paths
    ]

"""
READ PROTEIN FILES
"""
prot_paths = []
samples_prot = []

if protein_dir:
    prot_paths = glob.glob(os.path.join(protein_dir, '*.faa'))

    samples_prot = [
        os.path.splitext(os.path.basename(fp))[0]
        for fp in prot_paths
    ]


"""
DEBUG
"""
print(f"Genome samples: {samples_genome}")
print(f"Protein samples: {samples_prot}")


"""
PATTERNS
"""
PATTERN_LONG = "{sample}.fasta"
PATTERN_PROT = "{sample}.faa"


"""
MERGE SAMPLE SPACE
"""
samples_names = sorted(set(samples_genome) | set(samples_prot))

try:
    if config['args']['output'] is None:
        dir_out = os.path.join('sphae.out')
    else:
        dir_out = config['args']['output']
except KeyError:
    dir_out = os.path.join('sphae.out')

dir_annot = os.path.join(dir_out, 'PROCESSING', "genome-annotate")
dir_tree = os.path.join(dir_out, 'PROCESSING', "trees")
dir_final = os.path.join(dir_out, 'final-annotate')
dir_log = os.path.join(dir_out, 'logs')

dir_env = os.path.join(workflow.basedir, "envs")
dir_script = os.path.join(workflow.basedir, "scripts")


"""
MERGE SAMPLE SPACE
"""
samples_names = sorted(set(samples_genome) | set(samples_prot))

"""
LOG DIRECTORY
"""
dir = {'log': os.path.join(dir_out, 'logs')}

"""
LOG HELPER
"""
def copy_log_file():
    files = glob.glob(os.path.join(".snakemake", "log", "*.snakemake.log"))
    if not files:
        return None
    current_log = max(files, key=os.path.getmtime)
    target_log = os.path.join(dir['log'], os.path.basename(current_log))
    shutil.copy(current_log, target_log)

"""
ONSTART / SUCCESS / ERROR
"""
onstart:
    if os.path.isdir(dir["log"]):
        oldLogs = filter(re.compile(r'.*.log').match, os.listdir(dir["log"]))
        for logfile in oldLogs:
            os.unlink(os.path.join(dir["log"], logfile))

onsuccess:
    sys.stderr.write('\n\nSphae ran successfully!\n\n')

onerror:
    sys.stderr.write('\n\nSphae run failed\n\n')