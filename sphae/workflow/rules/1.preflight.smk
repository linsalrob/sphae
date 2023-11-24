import os
import attrmap as ap
import attrmap.utils as au
import glob
from metasnek import fastq_finder

"""
ONSTART/END/ERROR
"""
def copy_log_file():
    files = glob.glob(os.path.join(".snakemake", "log", "*.snakemake.log"))
    if not files:
        return None
    current_log = max(files, key=os.path.getmtime)
    #shell("cat " + current_log + " >> " + str(config.args.log))


onstart:
    """Cleanup old log files before starting"""
    if os.path.isdir(dir.log):
        oldLogs = filter(re.compile(r'.*.log').match, os.listdir(dir.log))
        for logfile in oldLogs:
            os.unlink(os.path.join(dir.log, logfile))

# onsuccess:
#     """Print a success message"""
#     sys.stderr.write('\n\nphage_genomes finished successfully!\n\n')
#     copy_log_file()

# onerror:
#     """Print an error message"""
#     sys.stderr.write('\n\nERROR: phage_genomes failed to finish.\n\n')
#     copy_log_file()


"""
DIRECTORIES
"""
dir = ap.AttrMap()
dir.out = config.args.output
dir.fastp = os.path.join(dir.out, 'PROCESSING' ,'results','fastp')
dir.nanopore = os.path.join(dir.out, 'PROCESSING','results','filtlong')
dir.assembly = os.path.join(dir.out, 'PROCESSING','assembly')
dir.megahit = os.path.join(dir.assembly, 'megahit')
dir.flye = os.path.join(dir.assembly, 'flye')
dir.genome = os.path.join(dir.out, 'PROCESSING','genome')
dir.cov = os.path.join(dir.out, 'PROCESSING','coverage')
dir.pharokka = os.path.join(dir.out, 'PROCESSING','pharokka')
dir.log = os.path.join(dir.out, 'logs')
dir.bench = os.path.join(dir.out, 'PROCESSING','bench')
dir.final = os.path.join(dir.out, 'RESULTS')

dir.env = os.path.join(workflow.basedir, "envs")
dir.script = os.path.join(workflow.basedir, "scripts")

# database dir
if config.args.db_dir is None:
    dir.db = os.path.join(workflow.basedir, 'databases')
else:
    dir.db = config.args.db_dir

# temp dir
if config.args.temp_dir is None:
    dir.temp = os.path.join(dir.out, "temp")
else:
    dir.temp = config.args.temp_dir


"""
PARSE SAMPLES

samples (dict):
    reads (dict):
        r1 (str): filepath
        r2 (str): filepath
    names (list): keys(samples["reads"])
"""

samples = ap.AttrMap()

samples.reads = fastq_finder.parse_samples_to_dictionary(config.args._input)
samples.names = list(au.get_keys(samples.reads))
samples = au.convert_state(samples, read_only=True)
fastq_finder.write_samples_tsv(samples.reads, os.path.join(dir.out, "samples.tsv"))


"""
Wildcard constraints
"""
wildcard_constraints:
    sample="[a-zA-Z0-9_-]+",
    host = ".{0}|\.hostRm|\.hostRm_s",
    subsample = ".{0}|\.subsampled"