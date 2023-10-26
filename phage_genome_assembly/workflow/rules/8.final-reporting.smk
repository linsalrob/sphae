"""
Summarizing the results to one directory
"""

rule summarize_paired:
    input:
        genomes= os.path.join(dir.genome, "{sample}-pr", "{sample}.fasta"),
        gbk=os.path.join(dir.pharokka, "{sample}-pr", "phynteny", "phynteny.gbk"),
        amr =os.path.join(dir.pharokka, "{sample}-pr", "top_hits_card.tsv"),
        vfdb=os.path.join(dir.pharokka, "{sample}-pr", "top_hits_vfdb.tsv"),
        plot=os.path.join(dir.pharokka, "{sample}-pr", "{sample}_pharokka_plot.png"),
    output:

    params:
        o=os.path.join(dir.out, 'RESULTS')