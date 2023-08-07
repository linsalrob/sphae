import attrmap as ap
import attrmap.utils as au


targets = ap.AttrMap()


targets.db = []

targets.db.append(os.path.join(dir.db, config.db.pfam_file))
targets.db.append(os.path.join(dir.db, 'pharokka_db', 'phrogs_db.index'))
targets.db.append(os.path.join(dir.db, 'checkv_db', 'checkv-db-v1.5', 'README.txt'))


targets.qc = []

if config.args.sequencing == 'paired':
    targets.qc.append(expand(os.path.join(dir.prinseq, "{sample}.paired.{r12}.fastq.gz"), sample=samples.names, r12=["R1","R2","S"]))
elif config.args.sequencing == 'longread':
    targets.qc.append(expand(os.path.join(dir.nanopore, "{sample}.single.fastq.gz"), sample=samples.names))


targets.assemble = []

if config.args.sequencing == 'paired':
    targets.assemble.append(expand(os.path.join(dir.megahit, "{sample}", "final.contigs.fa"),sample=samples.names))
    targets.assemble.append(expand(os.path.join(dir.megahit, "{sample}", "final.fastg"),sample=samples.names))
    targets.assemble.append(expand(os.path.join(dir.megahit, "{sample}", "results", "sample_coverage.tsv"),sample=samples.names))
    targets.assemble.append(expand(os.path.join(dir.megahit, "{sample}", "final.contigs_result_table.csv"),sample=samples.names))
    targets.assemble.append(expand(os.path.join(dir.megahit, "{sample}", "graph_seq_details_megahit.tsv"),sample=samples.names))
    targets.assemble.append(expand(os.path.join(dir.megahit, "{sample}", "checkv", "quality_summary.tsv"), sample=samples.names))
    targets.assemble.append(expand(os.path.join(dir.assembly, "{sample}-assembly-stats_megahit.csv"),sample=samples.names))
    targets.assemble.append(expand(os.path.join(dir.genome, "{sample}", "{sample}.fasta"), sample=samples.names))
elif config.args.sequencing == 'longread':
    targets.assemble.append(expand(os.path.join(dir.flye, "{sample}", "assembly{file}"),sample=samples.names, file=[".fasta", "_graph.gfa", "_info.txt"]))
    targets.assemble.append(expand(os.path.join(dir.flye, "{sample}", "consensus.fasta"),sample=samples.names))
    targets.assemble.append(expand(os.path.join(dir.flye, "{sample}", "results", "sample_coverage.tsv"),sample=samples.names))
    targets.assemble.append(expand(os.path.join(dir.flye,"{sample}","{file}"),sample=samples.names, file=["assembly_result_table.csv","graph_seq_details_flye.tsv"]))
    targets.assemble.append(expand(os.path.join(dir.flye, "{sample}", "checkv", "quality_summary.tsv"), sample=samples.names))
    targets.assemble.append(expand(os.path.join(dir.assembly,"{sample}-assembly-stats_flye.csv"),sample=samples.names))
    targets.assemble.append(expand(os.path.join(dir.genome, "{sample}", "{sample}.fasta"), sample=samples.names))


targets.annotate = []

targets.annotate.append(expand(os.path.join(dir.pharokka, "{sample}", "pharokka.gff"), sample=samples.names))

targets.coverage = []

if config.args.sequencing == 'paired':
    cov_dir = dir.megahit
    targets.coverage.append(expand(os.path.join(cov_dir, "{sample}", "results", "sample_coverage.tsv"), sample=samples.names))
    targets.coverage.append(expand(os.path.join(cov_dir, "{sample}", "temp","{sample}.bam"), sample=samples.names))
    targets.coverage.append(expand(os.path.join(cov_dir, "{sample}", "temp","{sample}.bam.bai"), sample=samples.names))
    targets.coverage.append(expand(os.path.join(cov_dir, "{sample}", "temp", "{sample}.gencov.tsv"), sample=samples.names))

elif config.args.sequencing == 'longread':
    cov_dir = dir.flye
    targets.coverage.append(expand(os.path.join(cov_dir, "{sample}", "results", "sample_coverage.tsv"), sample=samples.names))
    targets.coverage.append(expand(os.path.join(cov_dir, "{sample}", "temp","{sample}.bam"), sample=samples.names))
    targets.coverage.append(expand(os.path.join(cov_dir, "{sample}", "temp","{sample}.bam.bai"), sample=samples.names))
    targets.coverage.append(expand(os.path.join(cov_dir, "{sample}", "temp", "{sample}.gencov.tsv"), sample=samples.names))
