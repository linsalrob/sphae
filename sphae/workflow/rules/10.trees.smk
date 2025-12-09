rule get_terL:
    params:
        folder   = os.path.join(dir_annot, "phynteny-pr"),  # the directory
        interest = "terminase large subunit",
    output:
        os.path.join(dir_tree, "marker_proteins", "all_terL.faa"),
    localrule: True
    script:
        os.path.join(dir_script, "get_marker_ids.py")

rule get_portal:
    params:
        folder   = os.path.join(dir_annot, "phynteny-pr"),
        interest = "portal protein",
    output:
        os.path.join(dir_tree, "marker_proteins", "all_portal.faa"),
    localrule: True
    script:
        os.path.join(dir_script, "get_marker_ids.py")

rule msa_terL:
    input:
        os.path.join(dir_tree, "marker_proteins", "all_terL.faa")
    conda:
        os.path.join(dir_env, "trees.yaml")
    output:
        os.path.join(dir_tree, "marker_proteins", "all_terL.aln")
    shell:
        """
        if [[ -f {input} && -s {input} ]]; then
            mafft {input} > {output}
        else
            echo "Input {input} is empty or does not exist. Skipping alignment."
        fi
        touch {output}
        """

rule msa_portal:
    input:
        os.path.join(dir_tree, "marker_proteins", "all_portal.faa"),
    output:
        os.path.join(dir_tree, "marker_proteins", "all_portal.aln")
    conda:
        os.path.join(dir_env, "trees.yaml")
    shell:
        """
        if [[ -f {input} && -s {input} ]]; then
            mafft {input} > {output}
        else
            echo "Input {input} is empty or does not exist. Skipping alignment."
        fi
        touch {output}
        """

rule iqtree_terL:
    input:
        os.path.join(dir_tree, "marker_proteins", "all_terL.aln")
    output:
        os.path.join(dir_tree, "marker_proteins", "all_terL.nwk")
    params:
        prefix="terL"
    conda:
        os.path.join(dir_env, "trees.yaml")
    threads:
        config['resources']['bigjob']['cpu']
    shell:
        """
        if [[ -f {input} && -s {input} ]]; then
            cat {input}
            fasttree -nopr {input} > {output}
        else
            echo "Input {input} is empty or does not exist. Skipping tree generation."
        fi
        touch {output}
        """

rule iqtree_portal:
    input:
        os.path.join(dir_tree, "marker_proteins", "all_portal.aln")
    output:
        os.path.join(dir_tree, "marker_proteins", "all_portal.nwk")
    params:
        prefix="terL"
    conda:
        os.path.join(dir_env, "trees.yaml")
    threads:
        config['resources']['bigjob']['cpu']
    shell:
        """
        if [[ -f {input} && -s {input} ]]; then
            cat {input}
            fasttree -nopr {input} > {output}
        else
            echo "Input {input} is empty or does not exist. Skipping tree generation."
        fi
        touch {output}
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
        
        if [[ -f {input.terl} && -s {input.terl} ]]; then
            cp {input.terl} {output.terl}
        else
            echo "Input {input.terl} is empty or does not exist. Touching {output.terl}."
            touch {output.terl}
        fi

        if [[ -f {input.portal} && -s {input.portal} ]]; then
            cp {input.portal} {output.portal}
        else
            echo "Input {input.portal} is empty or does not exist. Touching {output.portal}."
            touch {output.portal}
        fi
        """
