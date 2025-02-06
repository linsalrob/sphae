import os, glob

SAMPLES = [os.path.splitext(os.path.basename(f))[0].rsplit('_R1', 1)[0] 
           for f in glob.glob(os.path.join(input_dir, '*.fasta*'))] if input_dir else []

#terminase large subunit
rule get_marker_proteins_terL:
    input:
        gbk=expand(os.path.join(dir_annot,"{sample}-phynteny", "phynteny.gbk"),sample=SAMPLES),
    params:
        interest="terminase large subunit",
        script=os.path.join(dir_script, 'get_marker_genes.py'),
    output:
        marker_proteins = os.path.join(dir_tree, "marker_proteins", "all_terL.faa")
    shell:
        """
        touch {output.marker_proteins}
        python3 {params.script} --genbank {input.gbk} \
            --protein "{params.interest}" \
            --output {output.marker_proteins}
        """

#portal protein
rule get_marker_proteins_portal:
    input:
        gbk=expand(os.path.join(dir_annot,"{sample}-phynteny", "phynteny.gbk"),sample=SAMPLES),
    params:
        interest="portal protein",
        script=os.path.join(dir_script, 'get_marker_genes.py'),
    output:
        marker_proteins = os.path.join(dir_tree, "marker_proteins", "all_portal.faa")
    shell:
        """
        touch {output.marker_proteins}
        python3 {params.script} --genbank {input.gbk} \
            --protein "{params.interest}" \
            --output {output.marker_proteins}
        """

rule msa_terL:
    input:
        marker_proteins = os.path.join(dir_tree, "marker_proteins", "all_terL.faa"),
    params:
        dirs= os.path.join(dir_tree, "marker_proteins"),
    output:
        os.path.join(dir_tree, "marker_proteins", "all_terL.aln")
    threads:
        config['resources']['bigjob']['cpu']
    container:
        "docker://biocontainers/mafft:v7.407-2-deb_cv1"
    shell:
        """
        mafft {input.marker_proteins} > {output} 
        """

rule msa_portal:
    input:
        marker_proteins = os.path.join(dir_tree, "marker_proteins", "all_portal.faa"),
    params:
        dirs= os.path.join(dir_tree, "marker_proteins"),
    output:
        os.path.join(dir_tree, "marker_proteins", "all_portal.aln")
    threads:
        config['resources']['bigjob']['cpu']
    container:
        "docker://staphb/mafft:7.526"
    shell:
        """
        mafft {input.marker_proteins} > {output} 
        """

rule iqtree_terL:
    input:
        os.path.join(dir_tree, "marker_proteins", "all_terL.aln")
    output:
        os.path.join(dir_tree, "marker_proteins", "all_terL.nwk")
    params:
        prefix="terL"
    container:
        "docker://biocontainers/fasttree:v2.1.10-2-deb_cv1"
    threads:
        config['resources']['bigjob']['cpu']
    shell:
        """
        fasttree -nopr {input} > {output} 
        """

rule iqtree_portal:
    input:
        os.path.join(dir_tree, "marker_proteins", "all_portal.aln")
    output:
        os.path.join(dir_tree, "marker_proteins", "all_portal.nwk")
    params:
        prefix="terL"
    container:
        "docker://biocontainers/fasttree:v2.1.10-2-deb_cv1"
    threads:
        config['resources']['bigjob']['cpu']
    shell:
        """
        fasttree -nopr {input} > {output} 
        """

rule_tree_output:
    input:
        terl=os.path.join(dir_final, "marker_proteins", "all_terL.nwk"),
        portal=os.path.join(dir_final, "marker_proteins", "all_portal.nwk")
    output:
        terl=os.path.join(dir_final, "trees", "all_terL.nwk"),
        portal=os.path.join(dir_final, "trees", "all_portal.nwk")
    shell:
        """
        cp {input.terl} {output.terl}
        cp {input.portal} {output.portal}
        """