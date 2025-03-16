"""
Running pharokka for annotation 
"""
rule pharokka_megahit:
    """Annotate genomes with Pharokka"""
    input:
        os.path.join(dir_annotate, "{sample}-pr-genomes", "{sample}_1.fasta"),
    params:
        inputdir=os.path.join(dir_annotate, "{sample}-pr-genomes"),
        output=os.path.join(dir_annotate, "pharokka-pr"),
        db = config['args']['pharokka_db'],
        genes= config['params']['gene-predict'],
    output:
        gbk=os.path.join(dir_annotate, "pharokka-pr", "{sample}_1_pharokka", "{sample}_1.gbk"),
        card=os.path.join(dir_annotate, "pharokka-pr", "{sample}_1_pharokka", "top_hits_card.tsv"),
        vfdb=os.path.join(dir_annotate, "pharokka-pr", "{sample}_1_pharokka", "top_hits_vfdb.tsv"),
        spacers=os.path.join(dir_annotate, "pharokka-pr", "{sample}_1_pharokka", "{sample}_1_minced_spacers.txt"),
        taxa=os.path.join(dir_annotate, "pharokka-pr", "{sample}_1_pharokka", "{sample}_1_top_hits_mash_inphared.tsv"),
        cdden=os.path.join(dir_annotate, "pharokka-pr", "{sample}_1_pharokka", "{sample}_1_length_gc_cds_density.tsv"),
        cds=os.path.join(dir_annotate, "pharokka-pr", "{sample}_1_pharokka", "{sample}_1_cds_functions.tsv")
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
            for f in {params.inputdir}/*; do 
                data="$(basename "$f" .fasta)"
                pharokka.py \
                    -i "$f" \
                    -o {params.output}/"$data"_pharokka \
                    -d {params.db} \
                    -g {params.genes} \
                    -t {threads} \
                    -f -p "$data"\
                    2> {log}
            done 
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

rule pharokka_flye:
    """Annotate genomes with Pharokka"""
    input:
        os.path.join(dir_annotate, "{sample}-sr-genomes", "{sample}_1.fasta")
    params:
        inputdir=os.path.join(dir_annotate, "{sample}-sr-genomes"),
        output=os.path.join(dir_annotate, "pharokka-sr"),
        db = config['args']['pharokka_db'],
        genes= config['params']['gene-predict'],
    output:
        gbk=os.path.join(dir_annotate, "pharokka-sr", "{sample}_1_pharokka", "{sample}_1.gbk"),
        card=os.path.join(dir_annotate, "pharokka-sr", "{sample}_1_pharokka", "top_hits_card.tsv"),
        vfdb=os.path.join(dir_annotate, "pharokka-sr", "{sample}_1_pharokka", "top_hits_vfdb.tsv"),
        spacers=os.path.join(dir_annotate, "pharokka-sr", "{sample}_1_pharokka", "{sample}_1_minced_spacers.txt"),
        taxa=os.path.join(dir_annotate, "pharokka-sr", "{sample}_1_pharokka", "{sample}_1_top_hits_mash_inphared.tsv"),
        cdden=os.path.join(dir_annotate, "pharokka-sr", "{sample}_1_pharokka", "{sample}_1_length_gc_cds_density.tsv"),
        cds=os.path.join(dir_annotate, "pharokka-sr", "{sample}_1_pharokka", "{sample}_1_cds_functions.tsv")
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
            for f in {params.inputdir}/*; do 
                data="$(basename "$f" .fasta)"
                pharokka.py \
                    -i "$f" \
                    -o {params.output}/"$data"_pharokka \
                    -d {params.db} \
                    -g {params.genes} \
                    -t {threads} \
                    -f -p "$data"\
                    2> {log}
            done
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