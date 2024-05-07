"""
Running Phold to improve the anntoation 
https://github.com/gbouras13/phold
"""
rule phold_run_paired:
    input:
        gbk=os.path.join(dir.pharokka, "{sample}-pr", "{sample}.gbk")
    params:
        predict=os.path.join(dir.pharokka,"{sample}-pr-predict"),
        o=os.path.join(dir.pharokka,"{sample}-pr-phold"),
        prefix="{sample}",
        db=os.path.join(dir.db, "phold")
    output:
        gbk=os.path.join(dir.pharokka,"{sample}-pr-phold","{sample}.gbk")
    threads: 
        config.resources.smalljob.cpu
    conda:
        os.path.join(dir.env, "phold.yaml")
    resources:
        mem_mb = config.resources.smalljob.mem,
        time = config.resources.smalljob.time
    log:
        os.path.join(dir.log, "phold.{sample}.log")
    benchmark:
        os.path.join(dir.bench,"phold_megahit_{sample}.txt")
    shell:
        """
        if [[ -s {input.gbk} ]] ; then
            phold predict -i {input.gbk} -o {params.predict} -p {params.prefix} -t {threads} --cpu -d {params.db} -f
            phold compare -i {input.gbk} --predictions_dir {params.predict} -p {params.prefix} -o {params.o} -t {threads} -d {params.db} -f
        else
            touch {output.gbk}
        fi
        """

rule phold_run_longreads:
    input:
        gbk=os.path.join(dir.pharokka, "{sample}-sr", "{sample}.gbk")
    params:
        predict=os.path.join(dir.pharokka,"{sample}-sr-predict"),
        o=os.path.join(dir.pharokka,"{sample}-sr-phold"),
        prefix="{sample}",
        db=os.path.join(dir.db, "phold")
    output:
        gbk=os.path.join(dir.pharokka,"{sample}-sr-phold","{sample}.gbk"),
    threads: 
        config.resources.smalljob.cpu
    conda:
        os.path.join(dir.env, "phold.yaml")
    resources:
        mem_mb = config.resources.smalljob.mem,
        time = config.resources.smalljob.time
    log:
        os.path.join(dir.log, "phold.{sample}.log")
    benchmark:
        os.path.join(dir.bench,"phold_flye_{sample}.txt")
    shell:
        """
        if [[ -s {input.gbk} ]] ; then
            phold predict -i {input.gbk} -o {params.predict} -p {params.prefix} -t {threads} --cpu -d {params.db} -f
            phold compare -i {input.gbk} --predictions_dir {params.predict} -p {params.prefix} -o {params.o} -t {threads} -d {params.db} -f
        else
            touch {output.gbk}
        fi
        """