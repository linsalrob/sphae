

targets ={'qc':[], 'assemble':[], 'annotate':[]}

targets['qc'] = []
if config['args']['sequencing'] == 'paired':
    for sample in samples_names:
        targets['qc'].append(expand(os.path.join(dir_fastp, "{sample}_subsampled_{r12}.fastq.gz"), sample=sample, r12=["R1", "R2"]))
elif config['args']['sequencing'] == 'longread':
    for sample in samples_names:
        targets['qc'].append(expand(os.path.join(dir_nanopore, "{sample}_filt.fastq.gz"), sample=sample))
       
if config['args']['sequencing'] == 'paired':
    targets['assemble'].append(expand(os.path.join(dir_megahit, "{sample}-pr", "final.contigs.fa"),sample=samples_names))
    targets['assemble'].append(expand(os.path.join(dir_megahit, "{sample}-pr", "final.gfa"),sample=samples_names))
    targets['assemble'].append(expand(os.path.join(dir_megahit, "{sample}-pr", "final.contigs_result_table.csv"),sample=samples_names))
    targets['assemble'].append(expand(os.path.join(dir_megahit, "{sample}-pr", "graph_seq_details_megahit.tsv"),sample=samples_names))
    targets['assemble'].append(expand(os.path.join(dir_megahit, "{sample}-pr", "checkv", "quality_summary.tsv"), sample=samples_names))
    targets['assemble'].append(expand(os.path.join(dir_assembly, "{sample}-assembly-stats_megahit.csv"),sample=samples_names))
    targets['assemble'].append(expand(os.path.join(dir_genome, "{sample}-pr", "{sample}.fasta"), sample=samples_names))
elif config['args']['sequencing'] == 'longread':
    targets['assemble'].append(expand(os.path.join(dir_flye, "{sample}-sr", "assembly.fasta"),sample=samples_names, file=[".fasta", "_graph.gfa", "_info.txt"]))
    targets['assemble'].append(expand(os.path.join(dir_flye, "{sample}-sr", "consensus.fasta"),sample=samples_names))
    targets['assemble'].append(expand(os.path.join(dir_flye,"{sample}-sr","{file}"),sample=samples_names, file=["assembly_result_table.csv","graph_seq_details_flye.tsv"]))
    targets['assemble'].append(expand(os.path.join(dir_flye, "{sample}-sr", "checkv", "quality_summary.tsv"), sample=samples_names))
    targets['assemble'].append(expand(os.path.join(dir_assembly,"{sample}-assembly-stats_flye.csv"),sample=samples_names))
    targets['assemble'].append(expand(os.path.join(dir_genome, "{sample}-sr", "{sample}.fasta"), sample=samples_names))

if config['args']['sequencing'] == 'paired':
    targets['annotate'].append(expand(os.path.join(dir_pharokka, "{sample}-pharokka", "{sample}.gbk"), sample=samples_names))
    targets['annotate'].append(expand(os.path.join(dir_pharokka, "{sample}-pharokka", "top_hits_card.tsv"),sample=samples_names))
    targets['annotate'].append(expand(os.path.join(dir_pharokka, "{sample}-pharokka", "top_hits_vfdb.tsv"),sample=samples_names))
    targets['annotate'].append(expand(os.path.join(dir_pharokka, "{sample}-pharokka", "{sample}_minced_spacers.txt"),sample=samples_names))
    targets['annotate'].append(expand(os.path.join(dir_pharokka, "{sample}-pharokka", "{sample}_top_hits_mash_inphared.tsv"),sample=samples_names))
    targets['annotate'].append(expand(os.path.join(dir_pharokka, "{sample}-pharokka", "{sample}_length_gc_cds_density.tsv"), sample=samples_names))
    targets['annotate'].append(expand(os.path.join(dir_pharokka, "{sample}-pharokka", "{sample}_cds_functions.tsv"), sample=samples_names))
    targets['annotate'].append(expand(os.path.join(dir_pharokka, "{sample}-pr-phold", "{sample}.gbk"),sample=samples_names))
    targets['annotate'].append(expand(os.path.join(dir_pharokka, "{sample}-pr-phold", "sub_db_tophits", "acr_cds_predictions.tsv"), sample=samples_names))
    targets['annotate'].append(expand(os.path.join(dir_pharokka, "{sample}-pr-phold", "sub_db_tophits", "card_cds_predictions.tsv"), sample=samples_names))
    targets['annotate'].append(expand(os.path.join(dir_pharokka, "{sample}-pr-phold", "sub_db_tophits", "defensefinder_cds_predictions.tsv"), sample=samples_names))
    targets['annotate'].append(expand(os.path.join(dir_pharokka, "{sample}-pr-phold", "sub_db_tophits", "vfdb_cds_predictions.tsv"), sample=samples_names))
    targets['annotate'].append(expand(os.path.join(dir_pharokka, "{sample}-pr-phynteny", "phynteny.gbk"), sample=samples_names))
    targets['annotate'].append(expand(os.path.join(dir_pharokka, "{sample}-pr-phynteny", "pharokka_plot.png"), sample=samples_names))
    targets['annotate'].append(expand(os.path.join(dir_final, "{sample}-pr", "{sample}_summary.txt"), sample=samples_names))
    targets['annotate'].append(expand(os.path.join(dir_final, "{sample}-pr", "{sample}_summary.functions"), sample=samples_names))
elif config['args']['sequencing']== 'longread':
    targets['annotate'].append(expand(os.path.join(dir_pharokka, "{sample}-pharokka", "{sample}.gbk"), sample=samples_names))
    targets['annotate'].append(expand(os.path.join(dir_pharokka, "{sample}-pharokka", "top_hits_card.tsv"),sample=samples_names))
    targets['annotate'].append(expand(os.path.join(dir_pharokka, "{sample}-pharokka", "top_hits_vfdb.tsv"),sample=samples_names))
    targets['annotate'].append(expand(os.path.join(dir_pharokka, "{sample}-pharokka", "{sample}_minced_spacers.txt"),sample=samples_names))
    targets['annotate'].append(expand(os.path.join(dir_pharokka, "{sample}-pharokka", "{sample}_top_hits_mash_inphared.tsv"),sample=samples_names))
    targets['annotate'].append(expand(os.path.join(dir_pharokka, "{sample}-pharokka", "{sample}_length_gc_cds_density.tsv"), sample=samples_names))
    targets['annotate'].append(expand(os.path.join(dir_pharokka, "{sample}-pharokka", "{sample}_cds_functions.tsv"), sample=samples_names))
    targets['annotate'].append(expand(os.path.join(dir_pharokka, "{sample}-sr-phold", "{sample}.gbk"),sample=samples_names))
    targets['annotate'].append(expand(os.path.join(dir_pharokka, "{sample}-sr-phold", "sub_db_tophits", "acr_cds_predictions.tsv"), sample=samples_names))
    targets['annotate'].append(expand(os.path.join(dir_pharokka, "{sample}-sr-phold", "sub_db_tophits", "card_cds_predictions.tsv"), sample=samples_names))
    targets['annotate'].append(expand(os.path.join(dir_pharokka, "{sample}-sr-phold", "sub_db_tophits", "defensefinder_cds_predictions.tsv"), sample=samples_names))
    targets['annotate'].append(expand(os.path.join(dir_pharokka, "{sample}-sr-phold", "sub_db_tophits", "vfdb_cds_predictions.tsv"), sample=samples_names))
    targets['annotate'].append(expand(os.path.join(dir_pharokka, "{sample}-sr-phynteny", "phynteny.gbk"), sample=samples_names))
    targets['annotate'].append(expand(os.path.join(dir_pharokka, "{sample}-sr-phynteny", "pharokka_plot.png"), sample=samples_names))
    targets['annotate'].append(expand(os.path.join(dir_final, "{sample}-sr", "{sample}_summary.txt"), sample=samples_names))
    targets['annotate'].append(expand(os.path.join(dir_final, "{sample}-sr", "{sample}_summary.txt"), sample=samples_names))
