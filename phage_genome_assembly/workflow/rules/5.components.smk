"""
Getting the components from the assembly graph
"""
    
rule components_spades:
    input:
        contigs = os.path.join(dir.spades, "{sample}-pr", "contigs.fasta"),
        path = os.path.join(dir.spades, "{sample}-pr", "contigs.paths"),
        graph = os.path.join(dir.spades, "{sample}-pr", "assembly_graph_after_simplification.gfa"),
    output:
        os.path.join(dir.spades, "{sample}-pr", "graph_seq_details_spades.tsv")
    params:
        o = os.path.join(dir.spades, "{sample}-pr"),
        assembler = 'spades'
    conda:
        os.path.join(dir.env, "graph.yaml")
    log:
        os.path.join(dir.log, "components_spades.{sample}.log")
    benchmark:
        os.path.join(dir.bench, "components_spades.{sample}.txt")
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
        out = os.path.join(dir.flye, "{sample}-sr", "graph_seq_details_flye.tsv"),
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