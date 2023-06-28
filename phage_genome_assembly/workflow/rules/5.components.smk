"""
Getting the components from the assembly graph
"""
    
rule components_megahit:
    input:
        contigs = os.path.join(dir.megahit, "{sample}", "final.contigs.fa"),
        graph = os.path.join(dir.megahit, "{sample}", "final.fastg"),
    output:
        os.path.join(dir.megahit, "{sample}", "graph_seq_details_megahit.tsv")
    params:
        o = os.path.join(dir.megahit, "{sample}"),
        assembler = 'megahit'
    conda:
        os.path.join(dir.env, "graph.yaml")
    log:
        os.path.join(dir.log, "components_megahit.{sample}.log")
    benchmark:
        os.path.join(dir.bench, "components_megahit.{sample}.txt")
    threads:
        config.resources.job.cpu
    resources:
        mem_mb=config.resources.job.mem,
        time=config.resources.job.time
    script:
        os.path.join(dir.script, 'components.py')


rule components_flye_nano:
    input:
        contigs = os.path.join(dir.flye, "{sample}", "assembly.fasta"),
        graph = os.path.join(dir.flye, "{sample}", "assembly_graph.gfa"),
        path = os.path.join(dir.flye, "{sample}", "assembly_info.txt")
    output:
        out = os.path.join(dir.flye, "{sample}", "graph_seq_details_flye.tsv"),
    params:
        o = os.path.join(dir.flye, "{sample}"),
        assembler = 'flye'
    log:
        os.path.join(dir.log, "components_flye_nano.{sample}.log")
    benchmark:
        os.path.join(dir.bench, "components_flye_nano.{sample}.txt")
    conda:
        os.path.join(dir.env, "graph.yaml")
    threads:
        config.resources.job.cpu
    resources:
        mem_mb=config.resources.job.mem,
        time=config.resources.job.time
    script:
        os.path.join(dir.script, 'components.py')