"""
Running pharokka for anntoation 
"""
rule pharokka_megahit:
    """Annotate genomes with Pharokka"""
    input:
        os.path.join(dir.genome, "{sample}-pr", "{sample}.fasta"),
    params:
        o=os.path.join(dir.pharokka, "{sample}-pr"),
        db=os.path.join(dir.db, "pharokka_db"),
        sp="{sample}"
    output:
        gff=os.path.join(dir.pharokka, "{sample}-pr", "{sample}.gbk"),
        plot=os.path.join(dir.pharokka, "{sample}-pr", "{sample}_pharokka_plot.png"),
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
            -f -p {params.sp}\
            2> {log}
        
        pharokka_plotter.py -i {input} -n {params.sp}_pharokka_plot -o {params.o} -p {params.sp} -f
        """



rule pharokka_flye:
    """Annotate genomes with Pharokka"""
    input:
        os.path.join(dir.genome, "{sample}-sr", "{sample}.fasta")
    params:
        o=os.path.join(dir.pharokka, "{sample}-sr"),
        db=os.path.join(dir.db, "pharokka_db"),
        sp="{sample}"
    output:
        gff=os.path.join(dir.pharokka, "{sample}-sr", "{sample}.gbk"),
        plot=os.path.join(dir.pharokka, "{sample}-sr", "{sample}_pharokka_plot.png"),
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
            -f -p {params.sp}\
            2> {log}

        pharokka_plotter.py -i {input} -n {params.sp}_pharokka_plot -o {params.o} -p {params.sp} -f
        """