"""

Generating a table with the graph compents for each contig - nanopore assemblies 
"""

rule components_flye_nano:
    input:
        contigs = os.path.join(ASSEMBLY, "{sample}-flye", "assembly.fasta"),
        graph = os.path.join(ASSEMBLY, "{sample}-flye", "assembly_graph.gfa"),
        path = os.path.join(ASSEMBLY, "{sample}-flye", "assembly_info.txt")
    output:
        out = os.path.join(ASSEMBLY, "{sample}-flye", "graph_seq_details_flye.tsv"),
    params:
        o = os.path.join(ASSEMBLY, "{sample}-flye"),
        assembler = 'flye'
    log:
        os.path.join(logs, "components_flye_{sample}.log")
    conda: "../envs/graph.yaml"
    resources:
        mem_mb=64000
    script:
        os.path.join('..', 'scripts', 'components.py')
