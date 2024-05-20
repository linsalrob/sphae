import os
import glob
import yaml
from metasnek import fastq_finder

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

dir_fastp = os.path.join(dir_out, 'PROCESSING' ,'fastp')
dir_nanopore = os.path.join(dir_out, 'PROCESSING','filtlong')
dir_assembly = os.path.join(dir_out, 'PROCESSING','assembly')
dir_megahit = os.path.join(dir_assembly, 'megahit')
dir_flye = os.path.join(dir_assembly, 'flye')
dir_genome = os.path.join(dir_out, 'PROCESSING','genome')
dir_pharokka = os.path.join(dir_out, 'PROCESSING','annotate')
dir_log = os.path.join(dir_out, 'logs')
dir_final = os.path.join(dir_out, 'RESULTS')

dir_env = os.path.join(workflow.basedir, "envs")
dir_script = os.path.join(workflow.basedir, "scripts")

# database dir
if config['args']['db_dir'] is None:
    dir_db = os.path.join(workflow.basedir, 'databases')
else:
    dir_db = config['args']['db_dir']

# temp dir
if config['args']['temp_dir'] is None:
    dir_temp = os.path.join(dir_out, "temp")
else:
    dir_temp = config['args']['temp_dir']

#reading the input files
input_dir = config['args']['input']

# List of file paths matching the pattern
if config['args']['sequencing'] == 'paired':
    file_paths = glob.glob(os.path.join(input_dir, '*_R1*.fastq*'))
    samples_names = [os.path.splitext(os.path.basename(file_path))[0].rsplit('_R1', 1)[0] for file_path in file_paths]
    extn = [os.path.splitext(os.path.basename(file_path))[0].rsplit('_R1', 1)[1] + os.path.splitext(os.path.basename(file_path))[1] for file_path in file_paths]
elif config['args']['sequencing'] == 'longread':
    file_paths = glob.glob(os.path.join(input_dir, '*.fastq*'))
    samples_names, extn = zip(*(os.path.splitext(os.path.basename(file_path)) if '.' in os.path.basename(file_path) else (os.path.basename(file_path), '') for file_path in file_paths))
    
print(f"Samples are {samples_names}")
print(f"Extensions are {extn}")

FQEXTN = extn[0]
PATTERN_R1 = '{sample}_R1' + FQEXTN
PATTERN_R2 = '{sample}_R2' + FQEXTN
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