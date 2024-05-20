import os
import glob
import yaml
import re
import shutil
import sys

"""
CONFIG FILE
"""
configfile: os.path.join(workflow.basedir, "..", "config", "config.yaml")

"""
DIRECTORIES
"""
dir = {}
#declaring output file
try:
    if config['args']['output'] is None:
        dir_out = os.path.join('sphae.out')
    else:
        dir_out = config['args']['output']
except KeyError:
    dir_out = os.path.join('sphae.out')

# database dir
if config['args']['db_dir'] is None:
    dir_db = os.path.join(workflow.basedir, 'databases')
else:
    dir_db = config['args']['db_dir']

dir_annot = os.path.join(dir_out,'PROCESSING')
dir_final = os.path.join(dir_out,'final-annotate')
dir_log = os.path.join(dir_out, 'logs')

dir_env = os.path.join(workflow.basedir, "envs")
dir_script = os.path.join(workflow.basedir, "scripts")

#reading the input files
#print (config)
try:
    input_dir = config['args']['genome']
except KeyError:
    input_dir = None 

#reading the sample names
if input_dir:
    # Reading the sample names
    file_paths = glob.glob(os.path.join(input_dir, '*.fasta*'))
    if file_paths:
        samples_names, extn = zip(*(os.path.splitext(os.path.basename(file_path)) if '.' in os.path.basename(file_path) else (os.path.basename(file_path), '') for file_path in file_paths))
    else:
        samples_names, extn = [], []
else:
    file_paths = []
    samples_names, extn = [], []
    
print(f"Samples are {samples_names}")
print(f"Extensions are {extn}")

FQEXTN = extn[0]
PATTERN_LONG='{sample}'+FQEXTN

"""ONSTART/END/ERROR
Tasks to perform at various stages the start and end of a run.
"""
def copy_log_file():
    files = glob.glob(os.path.join(".snakemake", "log", "*.snakemake.log"))
    if not files:
        return None
    current_log = max(files, key=os.path.getmtime)
    target_log = os.path.join(dir['log'], os.path.basename(current_log))
    shutil.copy(current_log, target_log)

dir = {'log': os.path.join(dir_out, 'logs')}
onstart:
    """Cleanup old log files before starting"""
    if os.path.isdir(dir["log"]):
        oldLogs = filter(re.compile(r'.*.log').match, os.listdir(dir["log"]))
        for logfile in oldLogs:
            os.unlink(os.path.join(dir["log"], logfile))
onsuccess:
    """Print a success message"""
    sys.stderr.write('\n\nSphaehost ran successfully!\n\n')
onerror:
    """Print an error message"""
    sys.stderr.write('\n\nSphaehost run failed\n\n')