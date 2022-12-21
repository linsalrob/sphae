"""
Add your preflight checks as pure Python code here.
e.g. Configure the run, declare directories, validate the input files etc.
This preflight check to confirm the database filepaths 
"""

"""CONFIGURATION"""
configfile: os.path.join(workflow.basedir, '../', 'config', 'config.yaml')
configfile: os.path.join(workflow.basedir, '../', 'config', 'databases.yaml')

"""Setting some variable"""
OUTDIR = config['output']
print(f"Output files will be saved to directory, {OUTDIR}\n")
logs = os.path.join(OUTDIR,'logs')

"""CHECK IF CUSTOM DATABASE DIRECTORY"""
if config['customDatabaseDirectory'] is None:
    databaseDir = os.path.join(workflow.basedir, 'databases')
else:
    databaseDir = config['customDatabaseDirectory']
DATABASES = config['customDatabaseDirectory']
print(f"Databases are being saved in, {DATABASES} \n")


"""Checking the phage contigs directory """
GENOMES = config['phage']

if (os.path.exists(OUTDIR)==True):
    GEN, EXTN = glob_wildcards(os.path.join(OUTDIR, GENOMES, "{sample}.{ext}"))
    if (len(GEN) == 0):
        sys.stderr.write(f"We did not find any files in {GEN}. Is this the right read dir?\n")
        sys.stderr.write(f"If the files are there, but running into an error, check filepaths\n")
        sys.exit(0)
    if len(set(EXTN)) != 1:
        sys.stderr.write("FATAL: You have more than one type of file extension\n\t")
        sys.stderr.write("\n\t".join(set(EXTN)))
        sys.stderr.write("\nWe don't know how to handle these\n")
        sys.exit(0)
    else:
        print(f"Found the phages, sample: {GEN}")
            
        PHAGEDIR = os.path.join(OUTDIR, GENOMES)

"""ONSTART/END/ERROR
Tasks to perform at various stages the start and end of a run.
"""
onsuccess:
    """Print a success message"""
    sys.stderr.write('\n\nAnnotation successfully completed!\n\n')

onerror:
    """Print an error message"""
    sys.stderr.write('\n\nAnnotation failed.\n\n')

