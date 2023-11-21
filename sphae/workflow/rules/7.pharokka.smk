"""
Running pharokka for anntoation 
"""
rule rename_contigs_megahit:
    input:
        fin=os.path.join(dir.genome, "{sample}-pr", "{sample}.fasta"),
    params:
        s ="{sample}"
    output:
        out=os.path.join(dir.genome, "{sample}-pr", "{sample}_genome.fasta"),
        csv=os.path.join(dir.genome, "{sample}-pr", "{sample}_temp.csv")
    localrule: True
    log:
        os.path.join(dir.log, "rename-contigs.{sample}.log")
    script:
        os.path.join(dir.script, 'rename_genomes.py')


rule pharokka_megahit:
    """Annotate genomes with Pharokka"""
    input:
        os.path.join(dir.genome, "{sample}-pr", "{sample}_genome.fasta"),
    params:
        o=os.path.join(dir.pharokka, "{sample}-pr"),
        db=os.path.join(dir.db, "pharokka_db"),
        sp="{sample}"
    output:
        gbk=os.path.join(dir.pharokka, "{sample}-pr", "{sample}.gbk"),
        plot=os.path.join(dir.pharokka, "{sample}-pr", "{sample}_pharokka_plot.png"),
        card=os.path.join(dir.pharokka, "{sample}-pr", "top_hits_card.tsv"),
        vfdb=os.path.join(dir.pharokka, "{sample}-pr", "top_hits_vfdb.tsv"),
        spacers=os.path.join(dir.pharokka, "{sample}-pr", "{sample}_minced_spacers.txt"),
        taxa=os.path.join(dir.pharokka, "{sample}-pr", "{sample}_top_hits_mash_inphared.tsv"),
        reorient=os.path.join(dir.pharokka, "{sample}-pr", "{sample}_dnaapler_reoriented.fasta"),
        cdden=os.path.join(dir.pharokka, "{sample}-pr", "{sample}_length_gc_cds_density.tsv")
    conda:
        os.path.join(dir.env, "pharokka.yaml")
    threads:
        config.resources.smalljob.cpu
    resources:
        mem_mb = config.resources.smalljob.mem,
        time = config.resources.smalljob.time
    log:
        os.path.join(dir.log, "pharokka.{sample}.log")
    shell:
        """
        pharokka.py \
            -i {input} \
            -o {params.o} \
            -d {params.db} \
            -t {threads} \
            --dnaapler \
            -f -p {params.sp}\
            2> {log}
        touch {output.reorient}
        pharokka_plotter.py -i {input} -n {params.sp}_pharokka_plot -o {params.o} -p {params.sp} -f
        """

rule rename_contigs_flye:
    input:
        fin=os.path.join(dir.genome, "{sample}-sr", "{sample}.fasta"),
    params:
        s ="{sample}"
    output:
        out=os.path.join(dir.genome, "{sample}-sr", "{sample}_genome.fasta"),
        csv=os.path.join(dir.genome, "{sample}-sr", "{sample}_temp.csv")
    localrule: True
    log:
        os.path.join(dir.log, "rename-contigs.{sample}.log")
    script:
        os.path.join(dir.script, 'rename_genomes.py')

rule pharokka_flye:
    """Annotate genomes with Pharokka"""
    input:
        os.path.join(dir.genome, "{sample}-sr", "{sample}_genome.fasta")
    params:
        o=os.path.join(dir.pharokka, "{sample}-sr"),
        db=os.path.join(dir.db, "pharokka_db"),
        sp="{sample}"
    output:
        gbk=os.path.join(dir.pharokka, "{sample}-sr", "{sample}.gbk"),
        plot=os.path.join(dir.pharokka, "{sample}-sr", "{sample}_pharokka_plot.png"),
        card=os.path.join(dir.pharokka, "{sample}-sr", "top_hits_card.tsv"),
        vfdb=os.path.join(dir.pharokka, "{sample}-sr", "top_hits_vfdb.tsv"),
        spacers=os.path.join(dir.pharokka, "{sample}-sr", "{sample}_minced_spacers.txt"),
        taxa=os.path.join(dir.pharokka, "{sample}-sr", "{sample}_top_hits_mash_inphared.tsv"),
        reorient=os.path.join(dir.pharokka, "{sample}-sr", "{sample}_dnaapler_reoriented.fasta"),
        cdden=os.path.join(dir.pharokka, "{sample}-sr", "{sample}_length_gc_cds_density.tsv")
    conda:
        os.path.join(dir.env, "pharokka.yaml")
    threads:
        config.resources.smalljob.cpu
    resources:
        mem_mb = config.resources.smalljob.mem,
        time = config.resources.smalljob.time
    log:
        os.path.join(dir.log, "pharokka.{sample}.log")
    shell:
        """
         pharokka.py \
            -i {input} \
            -o {params.o} \
            -d {params.db} \
            -t {threads} \
            --dnaapler \
            -f -p {params.sp}\
            2> {log}

        touch {output.reorient}

        pharokka_plotter.py -i {input} -n {params.sp}_pharokka_plot -o {params.o} -p {params.sp} -f
        """