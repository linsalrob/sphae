"""
Add your preflight checks as pure Python code here.
e.g. Configure the run, declare directories, validate the input files etc.
This preflight check to confirm the database filepaths 
"""

"""CONFIGURATION"""
configfile: os.path.join(workflow.basedir, '../', 'config', 'config.yaml')
configfile: os.path.join(workflow.basedir, '../', 'config', 'databases.yaml')

"""CHECK IF CUSTOM DATABASE DIRECTORY"""
if config['customDatabaseDirectory'] is None:
    databaseDir = os.path.join(workflow.basedir, 'databases')
else:
    databaseDir = config['customDatabaseDirectory']
DATABASES = config['customDatabaseDirectory']
print(f"Databases are being saved in, {DATABASES} \n")


"""ONSTART/END/ERROR
Tasks to perform at various stages the start and end of a run.
"""
onsuccess:
    """Print a success message"""
    sys.stderr.write('\n\nDatabases successfully installed!\n\n')

onerror:
    """Print an error message"""
    sys.stderr.write('\n\nERROR: Databases was not installed.\n\n')

