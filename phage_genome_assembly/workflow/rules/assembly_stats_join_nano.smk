"""

Joining results from coverm, viralverify and graph compents to one file - nanopore files
"""

rule join_assembly_stats_unicycler:
    input:
        coverm = os.path.join(ASSEMBLY, "{sample}-nanopore-unicycler", "{sample}-contigs.tsv"),
        viralverify = os.path.join(ASSEMBLY, "{sample}-viralverify-unicycler_nano", "assembly_result_table.csv"),
        comp = os.path.join(ASSEMBLY, "{sample}-nanopore-unicycler", "graph_seq_details_unicycler.tsv")
    output:
        tsv = os.path.join(ASSEMBLY, "{sample}-assembly-stats_unicycler.tsv")
    conda: "../envs/graph.yaml"
    shell:
        """
            if [[ -s {input.coverm} ]]; then
                python phage_genome_assembly/workflow/scripts/joining_stats.py -c {input.coverm} -v {input.viralverify} -g {input.comp} -o {output.tsv}
            fi
        """

rule join_assembly_stats_flye:
    input:
        coverm = os.path.join(ASSEMBLY, "{sample}-flye", "{sample}-contigs.tsv"),
        viralverify = os.path.join(ASSEMBLY, "{sample}-viralverify-flye", "assembly_result_table.csv"),
        comp = os.path.join(ASSEMBLY, "{sample}-flye", "graph_seq_details_flye.tsv")
    output:
        tsv = os.path.join(ASSEMBLY, "{sample}-assembly-stats_flye.tsv")
    conda: "../envs/graph.yaml"
    shell:
        """
            if [[ -s {input.coverm} ]]; then
                python phage_genome_assembly/workflow/scripts/joining_stats.py -c {input.coverm} -v {input.viralverify} -g {input.comp} -o {output.tsv}
            fi
        """