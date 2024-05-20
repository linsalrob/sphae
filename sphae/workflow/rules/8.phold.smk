"""
Running Phold to improve the anntoation 
https://github.com/gbouras13/phold
"""
rule phold_run_paired:
    input:
        gbk=os.path.join(dir_pharokka, "{sample}-pharokka", "{sample}.gbk")
    params:
        predict=os.path.join(dir_pharokka,"{sample}-pr-predict"),
        o=os.path.join(dir_pharokka,"{sample}-pr-phold"),
        prefix="{sample}",
        db=os.path.join(dir_db, "phold")
    output:
        gbk=os.path.join(dir_pharokka,"{sample}-pr-phold","{sample}.gbk"),
        acr=os.path.join(dir_pharokka,"{sample}-pr-phold","sub_db_tophits", "acr_cds_predictions.tsv"),
        card=os.path.join(dir_pharokka,"{sample}-pr-phold","sub_db_tophits", "card_cds_predictions.tsv"),
        defense=os.path.join(dir_pharokka,"{sample}-pr-phold","sub_db_tophits", "defensefinder_cds_predictions.tsv"),
        vfdb=os.path.join(dir_pharokka,"{sample}-pr-phold","sub_db_tophits", "vfdb_cds_predictions.tsv")
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
            phold predict -i {input.gbk} -o {params.predict} -p {params.prefix} -t {threads} --cpu -d {params.db} -f 2> {log}
            phold compare -i {input.gbk} --predictions_dir {params.predict} -p {params.prefix} -o {params.o} -t {threads} -d {params.db} -f 2> {log}
        else
            touch {output.gbk}
        fi
        """

rule phold_run_longreads:
    input:
        gbk=os.path.join(dir_pharokka, "{sample}-pharokka", "{sample}.gbk")
    params:
        predict=os.path.join(dir_pharokka,"{sample}-sr-predict"),
        o=os.path.join(dir_pharokka,"{sample}-sr-phold"),
        prefix="{sample}",
        db=os.path.join(dir_db, "phold")
    output:
        gbk=os.path.join(dir_pharokka,"{sample}-sr-phold","{sample}.gbk"),
        acr=os.path.join(dir_pharokka,"{sample}-sr-phold","sub_db_tophits", "acr_cds_predictions.tsv"),
        card=os.path.join(dir_pharokka,"{sample}-sr-phold","sub_db_tophits", "card_cds_predictions.tsv"),
        defense=os.path.join(dir_pharokka,"{sample}-sr-phold","sub_db_tophits", "defensefinder_cds_predictions.tsv"),
        vfdb=os.path.join(dir_pharokka,"{sample}-sr-phold","sub_db_tophits", "vfdb_cds_predictions.tsv")
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
            phold predict -i {input.gbk} -o {params.predict} -p {params.prefix} -t {threads} --cpu -d {params.db} -f 2> {log}
            phold compare -i {input.gbk} --predictions_dir {params.predict} -p {params.prefix} -o {params.o} -t {threads} -d {params.db} -f 2> {log}
        else
            touch {output.gbk}
        fi
        """
