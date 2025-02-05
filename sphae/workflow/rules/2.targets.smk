

targets ={'qc':[], 'assemble':[], 'annotate':[]}

targets['qc'] = []
if config['args']['sequencing'] == 'paired':
    for sample in samples_names:
        targets['qc'].append(expand(os.path.join(dir_fastp, "{sample}_subsampled_{r12}.fastq.gz"), sample=sample, r12=["R1", "R2"]))
        targets['qc'].append(expand(os.path.join(dir_fastp, "{sample}_fastp.txt"), sample=sample))
elif config['args']['sequencing'] == 'longread':
    for sample in samples_names:
        targets['qc'].append(expand(os.path.join(dir_nanopore, "{sample}_filt.fastq.gz"), sample=sample))
        targets['qc'].append(expand(os.path.join(dir_nanopore, "{sample}_filt.txt"), sample=sample))
       
if config['args']['sequencing'] == 'paired':
    targets['assemble'].append(expand(os.path.join(dir_megahit, "{sample}-pr", "final.contigs.fa"),sample=samples_names))
    targets['assemble'].append(expand(os.path.join(dir_megahit, "{sample}-pr", "final.gfa"),sample=samples_names))
    targets['assemble'].append(expand(os.path.join(dir_megahit, "{sample}-pr", "final.contigs_result_table.csv"),sample=samples_names))
    targets['assemble'].append(expand(os.path.join(dir_megahit, "{sample}-pr", "graph_seq_details_megahit.tsv"),sample=samples_names))
    targets['assemble'].append(expand(os.path.join(dir_megahit, "{sample}-pr", "checkv", "quality_summary.tsv"), sample=samples_names))
    targets['assemble'].append(expand(os.path.join(dir_assembly, "{sample}-assembly-stats_megahit.csv"),sample=samples_names))
    targets['assemble'].append(expand(os.path.join(dir_genome, "{sample}-pr", "{sample}.fasta"), sample=samples_names))
    targets['assemble'].append(expand(os.path.join(dir_annotate, "{sample}-pr-genomes", "{sample}_1.fasta"), sample=samples_names))
elif config['args']['sequencing'] == 'longread':
    targets['assemble'].append(expand(os.path.join(dir_flye, "{sample}-sr", "assembly.fasta"),sample=samples_names, file=[".fasta", "_graph.gfa", "_info.txt"]))
    targets['assemble'].append(expand(os.path.join(dir_flye, "{sample}-sr", "consensus.fasta"),sample=samples_names))#    targets['assemble'].append(expand(os.path.join(dir_flye,"{sample}-sr","{file}"),sample=samples_names, file=["assembly_result_table.csv","graph_seq_details_flye.tsv"]))
    targets['assemble'].append(expand(os.path.join(dir_flye, "{sample}-sr", "checkv", "quality_summary.tsv"), sample=samples_names))
    targets['assemble'].append(expand(os.path.join(dir_assembly,"{sample}-assembly-stats_flye.csv"),sample=samples_names))
    targets['assemble'].append(expand(os.path.join(dir_genome, "{sample}-sr", "{sample}.fasta"), sample=samples_names))
    targets['assemble'].append(expand(os.path.join(dir_annotate,  "{sample}-sr-genomes", "{sample}_1.fasta"), sample=samples_names))


if config['args']['sequencing'] == 'paired':
    targets['annotate'].append(expand(os.path.join(dir_annotate, "pharokka-pr", "{sample}_1_pharokka", "{sample}_1.gbk"), sample=samples_names))
    targets['annotate'].append(expand(os.path.join(dir_annotate, "pharokka-pr", "{sample}_1_pharokka", "top_hits_card.tsv"),sample=samples_names))
    targets['annotate'].append(expand(os.path.join(dir_annotate, "pharokka-pr", "{sample}_1_pharokka", "top_hits_vfdb.tsv"),sample=samples_names))
    targets['annotate'].append(expand(os.path.join(dir_annotate, "pharokka-pr", "{sample}_1_pharokka", "{sample}_1_minced_spacers.txt"),sample=samples_names))
    targets['annotate'].append(expand(os.path.join(dir_annotate, "pharokka-pr", "{sample}_1_pharokka", "{sample}_1_top_hits_mash_inphared.tsv"),sample=samples_names))
    targets['annotate'].append(expand(os.path.join(dir_annotate, "pharokka-pr", "{sample}_1_pharokka", "{sample}_1_length_gc_cds_density.tsv"), sample=samples_names))
    targets['annotate'].append(expand(os.path.join(dir_annotate, "pharokka-pr", "{sample}_1_pharokka", "{sample}_1_cds_functions.tsv"), sample=samples_names))
    targets['annotate'].append(expand(os.path.join(dir_annotate, "phold-pr", "{sample}_1_phold", "{sample}_1.gbk"),sample=samples_names))
    targets['annotate'].append(expand(os.path.join(dir_annotate, "phold-pr", "{sample}_1_phold", "sub_db_tophits", "acr_cds_predictions.tsv"), sample=samples_names))
    targets['annotate'].append(expand(os.path.join(dir_annotate, "phold-pr", "{sample}_1_phold", "sub_db_tophits", "card_cds_predictions.tsv"), sample=samples_names))
    targets['annotate'].append(expand(os.path.join(dir_annotate, "phold-pr", "{sample}_1_phold", "sub_db_tophits", "defensefinder_cds_predictions.tsv"), sample=samples_names))
    targets['annotate'].append(expand(os.path.join(dir_annotate, "phold-pr", "{sample}_1_phold", "sub_db_tophits", "vfdb_cds_predictions.tsv"), sample=samples_names))
    targets['annotate'].append(expand(os.path.join(dir_annotate, "phynteny-pr", "{sample}_1_phynteny", "phynteny.gbk"), sample=samples_names))
    targets['annotate'].append(expand(os.path.join(dir_annotate, "phynteny-pr", "{sample}_1_phynteny", "plots", "{sample}_1.png"), sample=samples_names))
    targets['annotate'].append(expand(os.path.join(dir_final, "{sample}-pr", "{sample}_summary.txt"), sample=samples_names))
    targets['annotate'].append(expand(os.path.join(dir_final, "{sample}-pr", "{sample}_1_summary.functions"), sample=samples_names))
elif config['args']['sequencing']== 'longread':
    targets['annotate'].append(expand(os.path.join(dir_annotate, "pharokka-sr", "{sample}_1_pharokka", "{sample}_1.gbk"), sample=samples_names))
    targets['annotate'].append(expand(os.path.join(dir_annotate, "pharokka-sr", "{sample}_1_pharokka", "top_hits_card.tsv"),sample=samples_names))
    targets['annotate'].append(expand(os.path.join(dir_annotate, "pharokka-sr", "{sample}_1_pharokka", "top_hits_vfdb.tsv"),sample=samples_names))
    targets['annotate'].append(expand(os.path.join(dir_annotate, "pharokka-sr", "{sample}_1_pharokka", "{sample}_1_minced_spacers.txt"),sample=samples_names))
    targets['annotate'].append(expand(os.path.join(dir_annotate, "pharokka-sr", "{sample}_1_pharokka", "{sample}_1_top_hits_mash_inphared.tsv"),sample=samples_names))
    targets['annotate'].append(expand(os.path.join(dir_annotate, "pharokka-sr", "{sample}_1_pharokka", "{sample}_1_length_gc_cds_density.tsv"), sample=samples_names))
    targets['annotate'].append(expand(os.path.join(dir_annotate, "pharokka-sr", "{sample}_1_pharokka", "{sample}_1_cds_functions.tsv"), sample=samples_names))
    targets['annotate'].append(expand(os.path.join(dir_annotate, "phold-sr", "{sample}_1_phold", "{sample}_1.gbk"),sample=samples_names))
    targets['annotate'].append(expand(os.path.join(dir_annotate, "phold-sr", "{sample}_1_phold", "sub_db_tophits", "acr_cds_predictions.tsv"), sample=samples_names))
    targets['annotate'].append(expand(os.path.join(dir_annotate, "phold-sr", "{sample}_1_phold", "sub_db_tophits", "card_cds_predictions.tsv"), sample=samples_names))
    targets['annotate'].append(expand(os.path.join(dir_annotate, "phold-sr", "{sample}_1_phold", "sub_db_tophits", "defensefinder_cds_predictions.tsv"), sample=samples_names))
    targets['annotate'].append(expand(os.path.join(dir_annotate, "phold-sr", "{sample}_1_phold", "sub_db_tophits", "vfdb_cds_predictions.tsv"), sample=samples_names))
    targets['annotate'].append(expand(os.path.join(dir_annotate, "phynteny-sr", "{sample}_1_phynteny", "phynteny.gbk"), sample=samples_names))
    targets['annotate'].append(expand(os.path.join(dir_annotate, "phynteny-sr", "{sample}_1_phynteny", "plots", "{sample}_1.png"), sample=samples_names))
    targets['annotate'].append(expand(os.path.join(dir_final, "{sample}-sr", "{sample}_summary.txt"), sample=samples_names))
    targets['annotate'].append(expand(os.path.join(dir_final, "{sample}-sr", "{sample}_1_summary.functions"), sample=samples_names))

