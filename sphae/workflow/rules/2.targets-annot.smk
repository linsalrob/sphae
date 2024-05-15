import attrmap as ap
import attrmap.utils as au

targets = ap.AttrMap()

targets.annotate = []
#pharokka plot
targets.annotate.append(expand(os.path.join(dir.pharokka, "{sample}_pharokka", "{sample}.gbk"), sample=samples.names))
targets.annotate.append(expand(os.path.join(dir.pharokka, "{sample}_pharokka", "{sample}_pharokka_plot.png"), sample=samples.names))
targets.annotate.append(expand(os.path.join(dir.pharokka, "{sample}_pharokka", "top_hits_card.tsv"),sample=samples.names))
targets.annotate.append(expand(os.path.join(dir.pharokka, "{sample}_pharokka", "top_hits_vfdb.tsv"),sample=samples.names))
targets.annotate.append(expand(os.path.join(dir.pharokka, "{sample}_pharokka", "{sample}_minced_spacers.txt"),sample=samples.names))
targets.annotate.append(expand(os.path.join(dir.pharokka, "{sample}_pharokka", "{sample}_top_hits_mash_inphared.tsv"),sample=samples.names))
targets.annotate.append(expand(os.path.join(dir.pharokka, "{sample}_pharokka", "{sample}_length_gc_cds_density.tsv"), sample=samples.names))
targets.annotate.append(expand(os.path.join(dir.pharokka, "{sample}_pharokka", "{sample}_cds_functions.tsv"), sample=samples.names))

#phold results 
targets.annotate.append(expand(os.path.join(dir.pharokka, "{sample}-phold", "{sample}.gbk"),sample=samples.names))
targets.annotate.append(expand(os.path.join(dir.pharokka,"{sample}-phold", "sub_db_tophits", "acr_cds_predictions.tsv"), sample=sample.names))
targets.annotate.append(expand(os.path.join(dir.pharokka,"{sample}-phold", "sub_db_tophits", "card_cds_predictions.tsv"), sample=sample.names))
targets.annotate.append(expand(os.path.join(dir.pharokka,"{sample}-phold", "sub_db_tophits", "defensefinder_cds_predictions.tsv"), sample=sample.names))
targets.annotate.append(expand(os.path.join(dir.pharokka,"{sample}-phold", "sub_db_tophits", "vfdb_cds_predictions.tsv"), sample=sample.names))

#phynteny
targets.annotate.append(expand(os.path.join(dir.pharokka, "{sample}_phynteny", "phynteny.gbk"), sample=samples.names))

#summary write up
targets.annotate.append(expand(os.path.join(dir.final,"{sample}", "{sample}_phynteny.gbk"), sample=sample_names))
targets.annotate.append(expand(os.path.join(dir.final, "{sample}", "{sample}_top_hits_card.tsv"), sample=sample_names))
targets.annotate.append(expand(os.path.join(dir.final, "{sample}", "{sample}_top_hits_vfdb.tsv"), sample=sample_names))
targets.annotate.append(expand(os.path.join(dir.final,"{sample}", "{sample}_minced_spacers.txt"), sample=sample_names))
targets.annotate.append(expand(os.path.join(dir.final, "{sample}", "{sample}_top_hits_mash_inphared.tsv"), sample=sample_names))