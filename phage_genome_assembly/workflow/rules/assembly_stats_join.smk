"""

Joining results from coverm, viralverify and graph compents to one file
"""

rule join_assembly_stats_spades:
    input:
        coverm = os.path.join(ASSEMBLY, "{sample}-spades", "{sample}-contigs.tsv"),
        viralverify = os.path.join(ASSEMBLY, "{sample}-viralverify-spades", "contigs_result_table.csv"),
        comp = os.path.join(ASSEMBLY, "{sample}-spades", "graph_seq_details_spades.tsv")
    output:
        tsv = os.path.join(ASSEMBLY, "{sample}-assembly-stats_spades.tsv")
    conda: "../envs/graph.yaml"
    shell:
        """
            if [[ -s {input.coverm} ]]; then
                python phage_genome_assembly/workflow/scripts/joining_stats.py -c {input.coverm} -v {input.viralverify} -g {input.comp} -o {output.tsv}
            fi
        """


rule join_assembly_stats_megahit:
    input:
        coverm = os.path.join(ASSEMBLY, "{sample}-megahit", "{sample}-contigs.tsv"),
        viralverify = os.path.join(ASSEMBLY, "{sample}-viralverify-megahit", "{sample}.contigs_result_table.csv"),
        comp = os.path.join(ASSEMBLY, "{sample}-megahit", "graph_seq_details_megahit.tsv")
    output:
        tsv = os.path.join(ASSEMBLY, "{sample}-assembly-stats_megahit.tsv")
    conda: "../envs/graph.yaml"
    shell:
        """
            if [[ -s {input.coverm} ]]; then
                python phage_genome_assembly/workflow/scripts/joining_stats.py -c {input.coverm} -v {input.viralverify} -g {input.comp} -o {output.tsv}
            fi        
        """