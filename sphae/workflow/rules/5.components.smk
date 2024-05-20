"""
Getting the components from the assembly graph
"""
    
rule components_megahit:
    input:
        contigs = os.path.join(dir_megahit, "{sample}-pr", "final.contigs.fa"),
        path = os.path.join(dir_megahit, "{sample}-pr", "final.fastg"),
        graph = os.path.join(dir_megahit, "{sample}-pr", "final.fastg"),
    output:
        os.path.join(dir_megahit, "{sample}-pr", "graph_seq_details_megahit.tsv")
    params:
        o = os.path.join(dir_megahit, "{sample}-pr"),
        assembler = 'megahit'
    conda:
        os.path.join(dir_env, "graph.yaml")
    log:
        os.path.join(dir_log, "components_megahit.{sample}.log")
    threads:
        config['resources']['smalljob']['cpu']
    resources:
        mem_mb=config['resources']['smalljob']['mem'],
        time=config['resources']['smalljob']['time']
    script:
        os.path.join(dir_script, 'components.py')


rule components_flye_nano:
    input:
        contigs = os.path.join(dir_flye, "{sample}-sr", "assembly.fasta"),
        graph = os.path.join(dir_flye, "{sample}-sr", "assembly_graph.gfa"),
        path = os.path.join(dir_flye, "{sample}-sr", "assembly_info.txt")
    output:
        os.path.join(dir_flye, "{sample}-sr", "graph_seq_details_flye.tsv"),
    params:
        o = os.path.join(dir_flye, "{sample}-sr"),
        assembler = 'flye'
    log:
        os.path.join(dir_log, "components_flye_nano.{sample}.log")
    conda:
        os.path.join(dir_env, "graph.yaml")
    threads:
        config['resources']['smalljob']['cpu']
    resources:
        mem_mb=config['resources']['smalljob']['mem'],
        time=config['resources']['smalljob']['time']
    script:
        os.path.join(dir_script, 'components.py')