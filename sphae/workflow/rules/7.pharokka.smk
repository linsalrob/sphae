"""
Running pharokka for anntoation 
"""
"""
Running pharokka for anntoation 
"""
rule rename_contigs_megahit:
    input:
        fin=os.path.join(dir_genome, "{sample}-pr", "{sample}.fasta"),
    params:
        s ="{sample}"
    output:
        out=os.path.join(dir_genome, "{sample}-pr", "{sample}_genome.fasta"),
        csv=os.path.join(dir_genome, "{sample}-pr", "{sample}_temp.csv")
    localrule: True
    log:
        os.path.join(dir_log, "rename-contigs.{sample}.log")
    script:
        os.path.join(dir_script, 'rename_genomes.py')

rule pharokka_megahit:
    """Annotate genomes with Pharokka"""
    input:
        os.path.join(dir_genome, "{sample}-pr", "{sample}_genome.fasta"),
    params:
        o=os.path.join(dir_pharokka, "{sample}-pharokka"),
        db=os.path.join(dir_db, "pharokka_db"),
        sp="{sample}"
    output:
        gbk=os.path.join(dir_pharokka, "{sample}-pharokka", "{sample}.gbk"),
        card=os.path.join(dir_pharokka, "{sample}-pharokka", "top_hits_card.tsv"),
        vfdb=os.path.join(dir_pharokka, "{sample}-pharokka", "top_hits_vfdb.tsv"),
        spacers=os.path.join(dir_pharokka, "{sample}-pharokka", "{sample}_minced_spacers.txt"),
        taxa=os.path.join(dir_pharokka, "{sample}-pharokka", "{sample}_top_hits_mash_inphared.tsv"),
        cdden=os.path.join(dir_pharokka, "{sample}-pharokka", "{sample}_length_gc_cds_density.tsv"),
        cds=os.path.join(dir_pharokka, "{sample}-pharokka", "{sample}_cds_functions.tsv")
    conda:
        os.path.join(dir_env, "pharokka.yaml")
    threads:
        config['resources']['smalljob']['cpu']
    resources:
        mem_mb = config['resources']['smalljob']['mem'],
        time = config['resources']['smalljob']['time']
    log:
        os.path.join(dir_log, "pharokka.{sample}.log")
    shell:
        """
        if [[ -s {input} ]] ; then
            pharokka.py \
                -i {input} \
                -o {params.o} \
                -d {params.db} \
                -t {threads} \
                -f -p {params.sp}\
                2> {log}
            touch {output.gbk}
            touch {output.card}
            touch {output.vfdb}
            touch {output.spacers}
            touch {output.taxa}
            touch {output.cdden}
            touch {output.cds}
        else
            touch {output.gbk}
            touch {output.card}
            touch {output.vfdb}
            touch {output.spacers}
            touch {output.taxa}
            touch {output.cdden}
            touch {output.cds}
        fi
        """

rule rename_contigs_flye:
    input:
        fin=os.path.join(dir_genome, "{sample}-sr", "{sample}.fasta"),
    params:
        s ="{sample}"
    output:
        out=os.path.join(dir_genome, "{sample}-sr", "{sample}_genome.fasta"),
        csv=os.path.join(dir_genome, "{sample}-sr", "{sample}_temp.csv")
    localrule: True
    log:
        os.path.join(dir_log, "rename-contigs.{sample}.log")
    script:
        os.path.join(dir_script, 'rename_genomes.py')

rule pharokka_flye:
    """Annotate genomes with Pharokka"""
    input:
        os.path.join(dir_genome, "{sample}-sr", "{sample}_genome.fasta")
    params:
        o=os.path.join(dir_pharokka, "{sample}-pharokka"),
        db=os.path.join(dir_db, "pharokka_db"),
        sp="{sample}"
    output:
        gbk=os.path.join(dir_pharokka, "{sample}-pharokka", "{sample}.gbk"),
        card=os.path.join(dir_pharokka, "{sample}-pharokka", "top_hits_card.tsv"),
        vfdb=os.path.join(dir_pharokka, "{sample}-pharokka", "top_hits_vfdb.tsv"),
        spacers=os.path.join(dir_pharokka, "{sample}-pharokka", "{sample}_minced_spacers.txt"),
        taxa=os.path.join(dir_pharokka, "{sample}-pharokka", "{sample}_top_hits_mash_inphared.tsv"),
        cdden=os.path.join(dir_pharokka, "{sample}-pharokka", "{sample}_length_gc_cds_density.tsv"),
        cds=os.path.join(dir_pharokka, "{sample}-pharokka", "{sample}_cds_functions.tsv")
    conda:
        os.path.join(dir_env, "pharokka.yaml")
    threads:
        config['resources']['smalljob']['cpu']
    resources:
        mem_mb = config['resources']['smalljob']['mem'],
        time = config['resources']['smalljob']['time']
    log:
        os.path.join(dir_log, "pharokka.{sample}.log")
    shell:
        """
        if [[ -s {input} ]] ; then
            pharokka.py \
                -i {input} \
                -o {params.o} \
                -d {params.db} \
                -t {threads} \
                -f -p {params.sp}\
                2> {log}

            touch {output.gbk}
            touch {output.card}
            touch {output.vfdb}
            touch {output.spacers}
            touch {output.taxa}
            touch {output.cdden}
            touch {output.cds}
        else
            touch {output.gbk}
            touch {output.card}
            touch {output.vfdb}
            touch {output.spacers}
            touch {output.taxa}
            touch {output.cdden}
            touch {output.cds}
        fi
        """