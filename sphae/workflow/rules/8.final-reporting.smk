"""
Summarizing the results to one directory
"""
rule summarize_paired:
    input:
        assembly=os.path.join(dir.megahit, "{sample}-pr", "log"),
        genome=os.path.join(dir.genome, "{sample}-pr", "{sample}.fasta"),
        gbk=os.path.join(dir.pharokka, "{sample}-pr", "phynteny", "phynteny.gbk"),
        table= os.path.join(dir.genome, "{sample}-pr", "{sample}-genome-candidates.csv"),
        amr =os.path.join(dir.pharokka, "{sample}-pr", "top_hits_card.tsv"),
        vfdb=os.path.join(dir.pharokka, "{sample}-pr", "top_hits_vfdb.tsv"),
        spacers=os.path.join(dir.pharokka, "{sample}-pr", "{sample}_minced_spacers.txt"),
        plot=os.path.join(dir.pharokka, "{sample}-pr", "{sample}_pharokka_plot.png"),
        ph_taxa =os.path.join(dir.pharokka, "{sample}-pr", "{sample}_top_hits_mash_inphared.tsv"),
        cdden=os.path.join(dir.pharokka, "{sample}-pr", "{sample}_length_gc_cds_density.tsv"),
    output:
        summary=os.path.join(dir.final, "{sample}-pr", "{sample}_summary.txt")
    params:
        genomes= os.path.join(dir.final, "{sample}-pr", "{sample}_genome.fasta"),
        gbks=os.path.join(dir.final, "{sample}-pr", "{sample}.gbk"),
        plots=os.path.join(dir.final, "{sample}-pr", "{sample}_pharokka_plot.png"),
        outdir=os.path.join(dir.final),
        sample="{sample}",
    localrule: True
    log: 
        os.path.join(dir.log, "final_summary.{sample}.log")
    script:
        os.path.join(dir.script, 'summary.py')

rule summarize_longread:
    input:
        assembly=os.path.join(dir.flye, "{sample}-sr", "flye.log"),
        genome=os.path.join(dir.genome, "{sample}-sr", "{sample}.fasta"),
        gbk=os.path.join(dir.pharokka, "{sample}-sr", "phynteny", "phynteny.gbk"),
        table= os.path.join(dir.genome, "{sample}-sr", "{sample}-genome-candidates.csv"),
        amr =os.path.join(dir.pharokka, "{sample}-sr", "top_hits_card.tsv"),
        vfdb=os.path.join(dir.pharokka, "{sample}-sr", "top_hits_vfdb.tsv"),
        plot=os.path.join(dir.pharokka, "{sample}-sr", "{sample}_pharokka_plot.png"),
        spacers=os.path.join(dir.pharokka, "{sample}-sr", "{sample}_minced_spacers.txt"),
        ph_taxa =os.path.join(dir.pharokka, "{sample}-sr", "{sample}_top_hits_mash_inphared.tsv"),
        cdden=os.path.join(dir.pharokka, "{sample}-sr", "{sample}_length_gc_cds_density.tsv"),
    output:
        summary=os.path.join(dir.final, "{sample}-sr", "{sample}_summary.txt")
    params:
        genomes= os.path.join(dir.final, "{sample}-sr", "{sample}_genome.fasta"),
        gbks=os.path.join(dir.final, "{sample}-sr", "{sample}.gbk"),
        plots=os.path.join(dir.final, "{sample}-sr", "{sample}_pharokka_plot.png"),
        sample="{sample}",
        outdir=os.path.join(dir.final)
    localrule: True
    log: 
        os.path.join(dir.log, "final_summary.{sample}.log")
    script:
        os.path.join(dir.script, 'summary.py')

