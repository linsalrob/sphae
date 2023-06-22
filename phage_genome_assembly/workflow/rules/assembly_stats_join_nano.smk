"""

Joining results from coverm, viralverify and graph compents to one file - nanopore files
"""

rule join_assembly_stats_flye:
    input:
        coverm = os.path.join(ASSEMBLY, "{sample}-flye", "{sample}-contigs.tsv"),
        viralverify = os.path.join(ASSEMBLY, "{sample}-viralverify-flye", "assembly_result_table.csv"),
        comp = os.path.join(ASSEMBLY, "{sample}-flye", "graph_seq_details_flye.tsv")
    output:
        csv = os.path.join(ASSEMBLY, "{sample}-assembly-stats_flye.csv")
    conda: "../envs/graph.yaml"
    script:
        os.path.join('..', 'scripts', 'joining_stats.py')