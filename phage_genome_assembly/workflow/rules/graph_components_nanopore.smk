"""

Generating a table with the graph compents for each contig - nanopore assemblies 
"""

rule components_unicycler_nano:
    input:
        contigs= os.path.join(ASSEMBLY, "{sample}-nanopore-unicycler", "assembly.fasta"),
        graph= os.path.join(ASSEMBLY, "{sample}-nanopore-unicycler", "assembly.gfa")
    output:
        os.path.join(ASSEMBLY, "{sample}-nanopore-unicycler", "graph_seq_details_unicycler.tsv")
    params:
        o = os.path.join(ASSEMBLY, "{sample}-nanopore-unicycler")
    log:
        os.path.join(logs, "components_unicycler_{sample}.log")
    conda: "../envs/graph.yaml"
    resources:
        mem_mb=64000
    shell:
        """
            if [[ -s {input.contigs} ]]; then
                python phage_genome_assembly/workflow/scripts/components.py -a unicycler -c {input.contigs} -g {input.graph} -o {params.o} 2> {log}
            fi
        """

rule components_flye_nano:
    input:
        contigs = os.path.join(ASSEMBLY, "{sample}-flye", "assembly.fasta"),
        graph = os.path.join(ASSEMBLY, "{sample}-flye", "assembly_graph.gfa"),
        path= os.path.join(ASSEMBLY, "{sample}-flye", "assembly_info.txt")
    output:
        os.path.join(ASSEMBLY, "{sample}-flye", "graph_seq_details_flye.tsv"),
    params:
        o = os.path.join(ASSEMBLY, "{sample}-flye")
    log:
        os.path.join(logs, "components_flye_{sample}.log")
    conda: "../envs/graph.yaml"
    resources:
        mem_mb=64000
    shell:
        """
            if [[ -s {input.contigs} ]]; then
                python phage_genome_assembly/workflow/scripts/components.py -a flye -c {input.contigs} -p {input.path} -g {input.graph} -o {params.o} 2> {log}
            fi
        """
