"""
Add your preflight checks as pure Python code here.
e.g. Configure the run, declare directories, validate the input files etc.
This preflight check to confirm the database filepaths 
"""

"""CONFIGURATION"""
configfile: os.path.join(workflow.basedir, '../', 'config', 'config.yaml')
configfile: os.path.join(workflow.basedir, '../', 'config', 'databases.yaml')

"""
Setting the directory variables
"""

READDIR = config['input']
OUTDIR = config['output']
print(f"Output files will be saved to directory, {OUTDIR}\n")
logs = os.path.join(OUTDIR,'logs')

if config ['sequencing'] == 'paired':
    QCDIR = os.path.join(OUTDIR, 'prinseq')
    print(f"Illumina QC is run using prinseq, and the output files are saved to, {QCDIR}\n")
    SAMPLES,EXTENSIONS, = glob_wildcards(os.path.join(READDIR, '{sample}_R1{extn}'))
    if len(SAMPLES) == 0:
        sys.stderr.write(f"We did not find any fastq files in {SAMPLES}. Is this the right read dir?\n")
        sys.stderr.write(f"If the files are there, but running into an error, check filepaths\n")
        sys.exit(0)
    if len(set(EXTENSIONS)) != 1:
        sys.stderr.write("FATAL: You have more than one type of file extension\n\t")
        sys.stderr.write("\n\t".join(set(EXTENSIONS)))
        sys.stderr.write("\nWe don't know how to handle these\n")
        sys.exit(0)

    FQEXTN = EXTENSIONS[0]
    PATTERN_R1 = '{sample}_R1' + FQEXTN
    PATTERN_R2 = '{sample}_R2' + FQEXTN

elif config['sequencing'] == 'longread':
    QCDIR = os.path.join(OUTDIR, 'filtlong')
    print(f"Nanopore fastq files run through QC using filtlong, the outputs are saved to, {QCDIR}\n")
    SAMPLES,EXTENSIONS, =glob_wildcards(os.path.join(READDIR, '{sample}.{extn}'))
    if len(SAMPLES) ==0:
        sys.stderr.write(f"We did not find any fastq files in {SAMPLES2}. Is this the right read dir?\n")
        sys.stderr.write(f"If the files are there, but running into an error, check filepaths\n")
        sys.exit(0)
    if len(set(EXTENSIONS)) != 1:
        sys.stderr.write("FATAL: You have more than one type of file extension\n\t")
        sys.stderr.write("\n\t".join(set(EXTENSIONS2)))
        sys.stderr.write("\nWe don't know how to handle these\n")
        sys.exit(0)

    FQEXTN = EXTENSIONS[0]
    PATTERN = '{sample}.'+FQEXTN

ASSEMBLY = os.path.join(OUTDIR, 'assembly')
print(f"Saving the assemblies here, {ASSEMBLY}\n")

############################################################################
#looking through phage contigs 
###########################################################################


PHAGECONTIGS = config['contigs']
print(f"Phage contigs assembled, {PHAGECONTIGS}\n")
if (os.path.exists(OUTDIR)==True):
    PHAGE, EXTN = glob_wildcards(os.path.join(OUTDIR, PHAGECONTIGS, "{sample}.{ext}"))
    if (len(PHAGE) == 0):
        sys.stderr.write(f"We did not find any files in {PHAGE}. Is this the right read dir?\n")
        sys.stderr.write(f"If the files are there, but running into an error, check filepaths\n")
        sys.exit(0)
    if len(set(EXTN)) != 1:
        sys.stderr.write("FATAL: You have more than one type of file extension\n\t")
        sys.stderr.write("\n\t".join(set(EXTN)))
        sys.stderr.write("\nWe don't know how to handle these\n")
        sys.exit(0)

    else:
        print(f"Found the phage contigs in phage_finalset, sample: {PHAGE}")
        if config ['sequencing'] == 'paired':
            QCDIR = os.path.join(OUTDIR, 'prinseq')
            CONTIGS=(set(SAMPLES).intersection(PHAGE))
        elif config['sequencing'] == 'longread':
            QCDIR = os.path.join(OUTDIR, 'filtlong')
            CONTIGS=(set(SAMPLES).intersection(PHAGE))
            
        GENOMEDIR = os.path.join(OUTDIR, PHAGECONTIGS)

"""ONSTART/END/ERROR
Tasks to perform at various stages the start and end of a run.
"""
onsuccess:
    """Print a success message"""
    sys.stderr.write('\n\nSuccess!\n\n')

onerror:
    """Print an error message"""
    sys.stderr.write('\n\nFailed\n\n')

