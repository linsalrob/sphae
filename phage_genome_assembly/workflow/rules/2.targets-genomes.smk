"""
Declare your targets here!
A separate file is ideal if you have lots of target files to create, or need some python logic to determine
the targets to declare. This example shows targets that are dependent on the input file type.
"""

genomes=[]

#recircularize the phage contigs to start with terminase subunit 
genomes.append(expand(os.path.join(OUTDIR, "recircular", "{sample}-terminase.tsv"), sample=PHAGE))
genomes.append(expand(os.path.join(OUTDIR, "recircular", "{sample}.fasta"), sample=PHAGE))
genomes.append(expand(os.path.join(OUTDIR, "recircular-rc", "{sample}.fasta"),sample=PHAGE))


if config['sequencing'] == 'paired':
    #coverage depths across the phage contigs
    genomes.append(expand(os.path.join(OUTDIR, "coverage", "{sample}-illuminaReads.tsv"), sample=CONTIGS))
    genomes.append(expand(os.path.join(OUTDIR, "coverage", "{sample}-illuminaReads-bam", "coverm-genome.{sample}_good_out_R1.fastq.bam"), sample=CONTIGS))
    genomes.append(expand(os.path.join(OUTDIR, "coverage", "{sample}-illuminaReads-bam", "coverm-genome.{sample}_good_out_R1.fastq.bam.bai"), sample=CONTIGS))
    genomes.append(expand(os.path.join(OUTDIR, "coverage", "{sample}-illuminaReads-bam", "{sample}-Ill-bedtools-genomecov.tsv"), sample=CONTIGS))
elif config['sequencing'] == 'longread':
    genomes.append(expand(os.path.join(OUTDIR, "coverage", "{sample}-NanoReads.tsv"), sample=CONTIGS))
    genomes.append(expand(os.path.join(OUTDIR, "coverage", "{sample}-NanoReads-bam", "coverm-genome.{sample}-filtlong.fastq.bam"), sample=CONTIGS))
    genomes.append(expand(os.path.join(OUTDIR, "coverage", "{sample}-illuminaReads-bam", "coverm-genome.{sample}_good_out_R1.fastq.bai"), sample=CONTIGS))
    genomes.append(expand(os.path.join(OUTDIR, "coverage", "{sample}-NanoReads-bam", "{sample}-Nano-bedtools-genomecov.tsv"), sample=CONTIGS))
