"""
Getting the components from the assembly graph
"""
rule components_spades:
    input:
        contigs = os.path.join(ASSEMBLY, "{sample}-spades", "contigs.fasta"),
        graph = os.path.join(ASSEMBLY, "{sample}-spades", "assembly_graph_with_scaffolds.gfa"),
        path = os.path.join(ASSEMBLY, "{sample}-spades", "contigs.paths"),
    output:
        os.path.join(ASSEMBLY, "{sample}-spades", "graph_seq_details_spades.tsv")
    params:
        o = os.path.join(ASSEMBLY, "{sample}-spades"),
        assembler = 'spades'
    conda: "../envs/graph.yaml"
    log:
        os.path.join(logs, "components_spades_{sample}.log")
    resources:
        mem_mb=64000
    script:
        os.path.join('..', 'scripts', 'components.py')

    
rule components_megahit:
    input:
        contigs = os.path.join(ASSEMBLY, "{sample}-megahit", "{sample}.contigs.fa"),
        graph = os.path.join(ASSEMBLY, "{sample}-megahit", "{sample}.fastg"),
        path =  os.path.join(ASSEMBLY, "{sample}-megahit", "{sample}.fastg") # don't actually need this for components.py with megahit but to ensure script runs
    output:
        os.path.join(ASSEMBLY, "{sample}-megahit", "graph_seq_details_megahit.tsv")
    params:
        o = os.path.join(ASSEMBLY, "{sample}-megahit"),
        assembler = 'megahit'
    conda: "../envs/graph.yaml"
    log:
        os.path.join(logs, "components_megahit_{sample}.log")
    resources:
        mem_mb=64000
    script:
        os.path.join('..', 'scripts', 'components.py')
