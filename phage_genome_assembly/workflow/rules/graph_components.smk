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
        o = os.path.join(ASSEMBLY, "{sample}-spades")
    conda: "../envs/graph.yaml"
    log:
        os.path.join(logs, "components_spades_{sample}.log")
    resources:
        mem_mb=64000
    shell:
        """
            if [[ -s {input.contigs} ]]; then
                python phage_genome_assembly/workflow/scripts/components.py -a spades -c {input.contigs} -g {input.graph} -p {input.path} -o {params.o} 2> {log}
            fi
        """
rule components_megahit:
    input:
        contigs = os.path.join(ASSEMBLY, "{sample}-megahit", "{sample}.contigs.fa"),
        graph = os.path.join(ASSEMBLY, "{sample}-megahit", "{sample}.fastg")
    output:
        os.path.join(ASSEMBLY, "{sample}-megahit", "graph_seq_details_megahit.tsv")
    params:
        o = os.path.join(ASSEMBLY, "{sample}-megahit")
    conda: "../envs/graph.yaml"
    log:
        os.path.join(logs, "components_megahit_{sample}.log")
    resources:
        mem_mb=64000
    shell:
        """
            if [[ -s {input.contigs} ]]; then
                python phage_genome_assembly/workflow/scripts/components.py -a megahit -c {input.contigs} -g {input.graph} -o {params.o} 2> {log}
            fi
        """
