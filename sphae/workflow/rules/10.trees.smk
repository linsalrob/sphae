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
        os.path.join(dir_tree, "marker_proteins", "all_terL.faa")
    container:
        "docker://nanozoo/mafft:7.508--c9dd46e"
    output:
        os.path.join(dir_tree, "marker_proteins", "all_terL.aln")
    shell:
        """
        mafft {input} > {output} 
        """

rule msa_portal:
    input:
        os.path.join(dir_tree, "marker_proteins", "all_portal.faa"),
    output:
        os.path.join(dir_tree, "marker_proteins", "all_portal.aln")
    container:
        "docker://biocontainers/mafft:v7.407-2-deb_cv1"
    shell:
        """
        mafft {input} > {output} 
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

rule tree_output:
    input:
        terl=os.path.join(dir_tree, "marker_proteins", "all_terL.nwk"),
        portal=os.path.join(dir_tree, "marker_proteins", "all_portal.nwk")
    params:
        dors=os.path.join(dir_final, "trees")
    output:
        terl=os.path.join(dir_final, "trees", "all_terL.nwk"),
        portal=os.path.join(dir_final, "trees", "all_portal.nwk")
    shell:
        """
	    mkdir -p {params.dors}
        mv {input.terl} {output.terl}
        mv {input.portal} {output.portal}
        """
