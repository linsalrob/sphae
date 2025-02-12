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
            mv {input.terl} {output.terl}
        else
            echo "Input {input.terl} is empty or does not exist. Touching {output.terl}."
            touch {output.terl}
        fi

        if [[ -f {input.portal} && -s {input.portal} ]]; then
            mv {input.portal} {output.portal}
        else
            echo "Input {input.portal} is empty or does not exist. Touching {output.portal}."
            touch {output.portal}
        fi
        """
