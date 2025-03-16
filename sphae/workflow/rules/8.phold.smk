"""
Running Phold to improve the anntoation 
https://github.com/gbouras13/phold
"""
rule phold_run_paired:
    input:
        gbk=os.path.join(dir_annotate, "pharokka-pr", "{sample}_1_pharokka", "{sample}_1.gbk")
    params:
        inputdir=os.path.join(dir_annotate, "{sample}-pr-genomes"),
        idir=os.path.join(dir_annotate, "pharokka-pr"),
        predict=os.path.join(dir_annotate, "predict-pr"),
        o=os.path.join(dir_annotate, "phold-pr"),
        db = config['args']['phold_db'],
    output:
        gbk=os.path.join(dir_annotate, "phold-pr", "{sample}_1_phold", "{sample}_1.gbk"),
        acr=os.path.join(dir_annotate, "phold-pr", "{sample}_1_phold", "sub_db_tophits", "acr_cds_predictions.tsv"),
        card=os.path.join(dir_annotate, "phold-pr", "{sample}_1_phold", "sub_db_tophits", "card_cds_predictions.tsv"),
        defense=os.path.join(dir_annotate, "phold-pr", "{sample}_1_phold", "sub_db_tophits", "defensefinder_cds_predictions.tsv"),
        vfdb=os.path.join(dir_annotate, "phold-pr", "{sample}_1_phold", "sub_db_tophits", "vfdb_cds_predictions.tsv")
    threads: 
        config['resources']['smalljob']['cpu']
    conda:
        os.path.join(dir_env, "phold.yaml")
    resources:
        mem_mb = config['resources']['smalljob']['mem'],
        time = config['resources']['smalljob']['time']
    log:
        os.path.join(dir_log, "phold.{sample}.log")
    shell:
        """
        if [[ -s {input.gbk} ]] ; then
            for f in {params.inputdir}/*; do 
                data="$(basename "$f" .fasta)"
                phold predict -i {params.idir}/"$data"_pharokka/"$data".gbk -o {params.predict}/"$data"_predict -p "$data" -t {threads} --cpu -d {params.db} -f 2> {log}
                phold compare -i {params.idir}/"$data"_pharokka/"$data".gbk --predictions_dir {params.predict}/"$data"_predict -p "$data" -o {params.o}/"$data"_phold -t {threads} -d {params.db} -f 2> {log}
            done
        else
            touch {output.gbk}
            touch {output.acr}
            touch {output.card}
            touch {output.defense}
            touch {output.vfdb}
        fi
        """

rule phold_run_longreads:
    input:
        gbk=os.path.join(dir_annotate, "pharokka-sr", "{sample}_1_pharokka", "{sample}_1.gbk")
    params:
        inputdir=os.path.join(dir_annotate, "{sample}-sr-genomes"),
        idir=os.path.join(dir_annotate, "pharokka-sr"),
        predict=os.path.join(dir_annotate, "predict-sr"),
        o=os.path.join(dir_annotate, "phold-sr"),
        db = config['args']['phold_db'],
    output:
        gbk=os.path.join(dir_annotate, "phold-sr", "{sample}_1_phold", "{sample}_1.gbk"),
        acr=os.path.join(dir_annotate, "phold-sr", "{sample}_1_phold", "sub_db_tophits", "acr_cds_predictions.tsv"),
        card=os.path.join(dir_annotate, "phold-sr", "{sample}_1_phold", "sub_db_tophits", "card_cds_predictions.tsv"),
        defense=os.path.join(dir_annotate, "phold-sr", "{sample}_1_phold", "sub_db_tophits", "defensefinder_cds_predictions.tsv"),
        vfdb=os.path.join(dir_annotate, "phold-sr", "{sample}_1_phold", "sub_db_tophits", "vfdb_cds_predictions.tsv")
    threads: 
        config['resources']['smalljob']['cpu']
    conda:
        os.path.join(dir_env, "phold.yaml")
    resources:
        mem_mb = config['resources']['smalljob']['mem'],
        time = config['resources']['smalljob']['time']
    log:
        os.path.join(dir_log, "phold.{sample}.log")
    shell:
        """
        if [[ -s {input.gbk} ]] ; then
            for f in {params.inputdir}/*; do 
                data="$(basename "$f" .fasta)"
                phold predict -i {params.idir}/"$data"_pharokka/"$data".gbk -o {params.predict}/"$data"_predict -p "$data" -t {threads} --cpu -d {params.db} -f 2> {log}
                phold compare -i {params.idir}/"$data"_pharokka/"$data".gbk --predictions_dir {params.predict}/"$data"_predict -p "$data" -o {params.o}/"$data"_phold -t {threads} -d {params.db} -f 2> {log}
            done
        else
            touch {output.gbk}
            touch {output.acr}
            touch {output.card}
            touch {output.defense}
            touch {output.vfdb}
        fi
        """
