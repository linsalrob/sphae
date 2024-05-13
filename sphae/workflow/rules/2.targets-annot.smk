import attrmap as ap
import attrmap.utils as au

targets = ap.AttrMap()

targets.annotate = []
targets.annotate.append(expand(os.path.join(dir.pharokka, "{sample}", "{sample}.gbk"), sample=samples.names))
targets.annotate.append(expand(os.path.join(dir.pharokka, "{sample}", "{sample}_pharokka_plot.png"), sample=samples.names))
targets.annotate.append(expand(os.path.join(dir.pharokka, "{sample}", "phynteny", "phynteny.gbk"), sample=samples.names))
targets.annotate.append(expand(os.path.join(dir.pharokka, "{sample}", "top_hits_card.tsv"),sample=samples.names))
targets.annotate.append(expand(os.path.join(dir.pharokka, "{sample}", "top_hits_vfdb.tsv"),sample=samples.names))
targets.annotate.append(expand(os.path.join(dir.pharokka, "{sample}", "{sample}_minced_spacers.txt"),sample=samples.names))
targets.annotate.append(expand(os.path.join(dir.pharokka, "{sample}-phold", "{sample}.gbk"),sample=samples.names))
targets.annotate.append(expand(os.path.join(dir.pharokka, "{sample}", "phynteny", "phynteny.gbk"), sample=samples.names))
targets.annotate.append(expand(os.path.join(dir.pharokka, "{sample}", "{sample}_top_hits_mash_inphared.tsv"),sample=samples.names))
targets.annotate.append(expand(os.path.join(dir.pharokka, "{sample}", "{sample}_length_gc_cds_density.tsv"), sample=samples.names))
targets.annotate.append(expand(os.path.join(dir.pharokka, "{sample}", "{sample}_cds_functions.tsv"), sample=samples.names))
