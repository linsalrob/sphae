import attrmap as ap
# import attrmap.utils as au

targets = ap.AttrMap()

targets.qc = []
if config.args.sequencing == 'paired':
    for sample in samples.names:
        targets.qc.append(expand(os.path.join(dir.fastp, "{sample}_{r12}.subsampled.fastq.gz"), sample=sample, r12=["R1", "R2", "RS"]))
elif config.args.sequencing == 'longread':
    for sample in samples.names:
        targets.qc.append(expand(os.path.join(dir.nanopore, "{sample}_S.subsampled.fastq.gz"), sample=sample))
       
targets.assemble = []

if config.args.sequencing == 'paired':
    targets.assemble.append(expand(os.path.join(dir.megahit, "{sample}-pr", "final.contigs.fa"),sample=samples.names))
    targets.assemble.append(expand(os.path.join(dir.megahit, "{sample}-pr", "final.gfa"),sample=samples.names))
    targets.assemble.append(expand(os.path.join(dir.megahit, "{sample}-pr", "final.contigs_result_table.csv"),sample=samples.names))
    targets.assemble.append(expand(os.path.join(dir.megahit, "{sample}-pr", "graph_seq_details_megahit.tsv"),sample=samples.names))
    targets.assemble.append(expand(os.path.join(dir.megahit, "{sample}-pr", "checkv", "quality_summary.tsv"), sample=samples.names))
    targets.assemble.append(expand(os.path.join(dir.assembly, "{sample}-assembly-stats_megahit.csv"),sample=samples.names))
    targets.assemble.append(expand(os.path.join(dir.genome, "{sample}-pr", "{sample}.fasta"), sample=samples.names))
elif config.args.sequencing == 'longread':
    targets.assemble.append(expand(os.path.join(dir.flye, "{sample}-sr", "assembly.fasta"),sample=samples.names, file=[".fasta", "_graph.gfa", "_info.txt"]))
    targets.assemble.append(expand(os.path.join(dir.flye, "{sample}-sr", "consensus.fasta"),sample=samples.names))
    targets.assemble.append(expand(os.path.join(dir.flye,"{sample}-sr","{file}"),sample=samples.names, file=["assembly_result_table.csv","graph_seq_details_flye.tsv"]))
    targets.assemble.append(expand(os.path.join(dir.flye, "{sample}-sr", "checkv", "quality_summary.tsv"), sample=samples.names))
    targets.assemble.append(expand(os.path.join(dir.assembly,"{sample}-assembly-stats_flye.csv"),sample=samples.names))
    targets.assemble.append(expand(os.path.join(dir.genome, "{sample}-sr", "{sample}.fasta"), sample=samples.names))


targets.annotate = []
if config.args.sequencing == 'paired':
    targets.annotate.append(expand(os.path.join(dir.pharokka, "{sample}-pr", "{sample}.gbk"), sample=samples.names))
    targets.annotate.append(expand(os.path.join(dir.pharokka, "{sample}-pr", "{sample}_pharokka_plot.png"), sample=samples.names))
    targets.annotate.append(expand(os.path.join(dir.pharokka, "{sample}-pr", "top_hits_card.tsv"),sample=samples.names))
    targets.annotate.append(expand(os.path.join(dir.pharokka, "{sample}-pr", "top_hits_vfdb.tsv"),sample=samples.names))
    targets.annotate.append(expand(os.path.join(dir.pharokka, "{sample}-pr", "{sample}_minced_spacers.txt"),sample=samples.names))
    targets.annotate.append(expand(os.path.join(dir.pharokka, "{sample}-pr", "phynteny", "phynteny.gbk"), sample=samples.names))
    targets.annotate.append(expand(os.path.join(dir.pharokka, "{sample}-pr", "{sample}_top_hits_mash_inphared.tsv"),sample=samples.names))
    targets.annotate.append(expand(os.path.join(dir.final, "{sample}-pr", "{sample}_summary.txt"), sample=samples.names))
    targets.annotate.append(expand(os.path.join(dir.pharokka, "{sample}-pr", "{sample}_length_gc_cds_density.tsv"), sample=samples.names))
elif config.args.sequencing == 'longread':
    targets.annotate.append(expand(os.path.join(dir.pharokka, "{sample}-sr", "{sample}.gbk"), sample=samples.names))
    targets.annotate.append(expand(os.path.join(dir.pharokka, "{sample}-sr", "{sample}_pharokka_plot.png"), sample=samples.names))
    targets.annotate.append(expand(os.path.join(dir.pharokka, "{sample}-sr", "phynteny", "phynteny.gbk"), sample=samples.names))
    targets.annotate.append(expand(os.path.join(dir.pharokka, "{sample}-sr", "top_hits_card.tsv"),sample=samples.names))
    targets.annotate.append(expand(os.path.join(dir.pharokka, "{sample}-sr", "top_hits_vfdb.tsv"),sample=samples.names))
    targets.annotate.append(expand(os.path.join(dir.pharokka, "{sample}-sr", "{sample}_minced_spacers.txt"),sample=samples.names))
    targets.annotate.append(expand(os.path.join(dir.pharokka, "{sample}-sr", "phynteny", "phynteny.gbk"), sample=samples.names))
    targets.annotate.append(expand(os.path.join(dir.pharokka, "{sample}-sr", "{sample}_top_hits_mash_inphared.tsv"),sample=samples.names))
    targets.annotate.append(expand(os.path.join(dir.final, "{sample}-sr", "{sample}_summary.txt"), sample=samples.names))
    targets.annotate.append(expand(os.path.join(dir.pharokka, "{sample}-sr", "{sample}_length_gc_cds_density.tsv"), sample=samples.names))
