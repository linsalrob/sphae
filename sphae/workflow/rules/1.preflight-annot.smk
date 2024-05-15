import os
import attrmap as ap
import attrmap.utils as au
import glob

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

"""
DIRECTORIES
"""
dir = ap.AttrMap()
dir.out = config.args.output
dir.pharokka = os.path.join(dir.out,'PROCESSING')
dir.final = os.path.join(dir.out,'annotation')
dir.log = os.path.join(dir.out, 'logs')
dir.bench = os.path.join(dir.out, 'benchmarks')

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

#reading the sample names
samples = ap.AttrMap()
samples.reads = config.args._input
file_names = os.listdir(samples.reads)
sample_names = [os.path.splitext(file_name)[0] for file_name in file_names if file_name.endswith('.fasta')]
samples.names = sample_names
samples = au.convert_state(samples, read_only=True)
#print (samples)