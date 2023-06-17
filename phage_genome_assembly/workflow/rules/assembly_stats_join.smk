"""

Joining results from coverm, viralverify and graph compents to one file
"""
rule join_assembly_stats_megahit:
    input:
        coverm = os.path.join(ASSEMBLY, "{sample}-megahit", "{sample}-contigs.tsv"),
        viralverify = os.path.join(ASSEMBLY, "{sample}-viralverify-megahit", "{sample}.contigs_result_table.csv"),
        comp = os.path.join(ASSEMBLY, "{sample}-megahit", "graph_seq_details_megahit.tsv")
    output:
        csv = os.path.join(ASSEMBLY, "{sample}-assembly-stats_megahit.csv")
    conda: "../envs/graph.yaml"
    script:
        os.path.join('..', 'scripts', 'joining_stats.py')