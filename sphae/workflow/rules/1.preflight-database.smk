import attrmap as ap
#import attrmap.utils as au
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

# onsuccess:
#     """Print a success message"""
#     sys.stderr.write('\n\nSpae finished successfully!\n\n')
#     copy_log_file()

# onerror:
#     """Print an error message"""
#     sys.stderr.write('\n\nERROR: Spae failed to finish.\n\n')
#     copy_log_file()

"""
DIRECTORIES
"""
dir = ap.AttrMap()
dir.env = os.path.join(workflow.basedir, "envs")

# database dir
if config.args.db_dir is None:
    databaseDir = os.path.join(workflow.basedir, 'databases')
else:
    databaseDir = config.args.db_dir
print(f"Databases are being saved in, {databaseDir} \n")
