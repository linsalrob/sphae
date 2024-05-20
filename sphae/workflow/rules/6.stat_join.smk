rule join_assembly_stats_megahit:
    input:
        viralverify = os.path.join(dir_megahit, "{sample}-pr", "final.contigs_result_table.csv"),
        comp = os.path.join(dir_megahit, "{sample}-pr", "graph_seq_details_megahit.tsv"),
        checkv = os.path.join(dir_megahit, "{sample}-pr", "checkv", "quality_summary.tsv")
    output:
        csv = os.path.join(dir_assembly, "{sample}-assembly-stats_megahit.csv")
    conda:
        os.path.join(dir_env, "graph.yaml")
    script:
        os.path.join(dir_script, 'joining_stats.py')


rule join_assembly_stats_flye:
    input:
        viralverify = os.path.join(dir_flye, "{sample}-sr", "assembly_result_table.csv"),
        comp = os.path.join(dir_flye, "{sample}-sr", "graph_seq_details_flye.tsv"),
        checkv = os.path.join(dir_flye, "{sample}-sr", "checkv", "quality_summary.tsv")
    output:
        csv = os.path.join(dir_assembly, "{sample}-assembly-stats_flye.csv")
    conda:
        os.path.join(dir_env, "graph.yaml")
    script:
        os.path.join(dir_script, 'joining_stats.py')
