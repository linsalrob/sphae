"""
Getting the components from the assembly graph
"""
    
rule components_megahit:
    input:
        contigs = os.path.join(dir.megahit, "{sample}-pr", "final.contigs.fa"),
        path = os.path.join(dir.megahit, "{sample}-pr", "final.fastg"),
        graph = os.path.join(dir.megahit, "{sample}-pr", "final.fastg"),
    output:
        os.path.join(dir.megahit, "{sample}-pr", "graph_seq_details_megahit.tsv")
    params:
        o = os.path.join(dir.megahit, "{sample}-pr"),
        assembler = 'megahit'
    conda:
        os.path.join(dir.env, "graph.yaml")
    log:
        os.path.join(dir.log, "components_megahit.{sample}.log")
    benchmark:
        os.path.join(dir.bench, "components_megahit.{sample}.txt")
    threads:
        config.resources.smalljob.cpu
    resources:
        mem_mb=config.resources.smalljob.mem,
        time=config.resources.smalljob.time
    script:
        os.path.join(dir.script, 'components.py')


rule components_flye_nano:
    input:
        contigs = os.path.join(dir.flye, "{sample}-sr", "assembly.fasta"),
        graph = os.path.join(dir.flye, "{sample}-sr", "assembly_graph.gfa"),
        path = os.path.join(dir.flye, "{sample}-sr", "assembly_info.txt")
    output:
        os.path.join(dir.flye, "{sample}-sr", "graph_seq_details_flye.tsv"),
    params:
        o = os.path.join(dir.flye, "{sample}-sr"),
        assembler = 'flye'
    log:
        os.path.join(dir.log, "components_flye_nano.{sample}.log")
    benchmark:
        os.path.join(dir.bench, "components_flye_nano.{sample}.txt")
    conda:
        os.path.join(dir.env, "graph.yaml")
    threads:
        config.resources.smalljob.cpu
    resources:
        mem_mb=config.resources.smalljob.mem,
        time=config.resources.smalljob.time
    script:
        os.path.join(dir.script, 'components.py')