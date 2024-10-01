"""
Summarizing the results to one directory
"""
rule summarize_annotations_paired:
    input: 
        pharokka=os.path.join(dir_annotate, "pharokka-pr", "{sample}_1_pharokka", "{sample}_1.gbk"),
        phold=os.path.join(dir_annotate, "phold-pr", "{sample}_1_phold", "{sample}_1.gbk"),
        pkl=os.path.join(dir_annotate, "phynteny-pr", "{sample}_1_phynteny", "phynteny.gbk")
    params:
        inputdir=os.path.join(dir_annotate, "{sample}-pr-genomes"),
        pharokka=os.path.join(dir_annotate, "pharokka-pr"),
        phold=os.path.join(dir_annotate, "phold-pr"),
        pkl=os.path.join(dir_annotate, "phynteny-pr")
    output:
        pharokka_func=os.path.join(dir_annotate, "pharokka-pr", "{sample}_1_pharokka", "{sample}_1_pharokka.functions"),
        phold_func=os.path.join(dir_annotate, "phold-pr", "{sample}_1_phold", "{sample}_1_phold.functions"),
        pkl_func=os.path.join(dir_annotate, "phynteny-pr", "{sample}_1_phynteny", "phynteny.functions"),
    resources:
        mem =config['resources']['smalljob']['mem'],
        time = config['resources']['smalljob']['time']
    conda:
        os.path.join(dir_env, "pharokka.yaml")
    shell:
        """
        if [[ -s {input.pharokka} ]] ; then
            for f in {params.inputdir}/*; do 
                data="$(basename "$f" .fasta)"
                genbank_to -g {params.pharokka}/"$data"_pharokka/"$data".gbk -f {params.pharokka}/"$data"_pharokka/"$data"_pharokka.functions
            done
        else
            touch {output.pharokka_func}
        fi

        if [[ -s {input.phold} ]] ; then
            for f in {params.inputdir}/*; do 
                data="$(basename "$f" .fasta)"
                genbank_to -g {params.phold}/"$data"_phold/"$data".gbk -f {params.phold}/"$data"_phold/"$data"_phold.functions
            done
        else
            touch {output.phold_func}
        fi

        if [[ -s {input.pkl} ]] ; then
            for f in {params.inputdir}/*; do 
                data="$(basename "$f" .fasta)"
                genbank_to -g {params.pkl}/"$data"_phynteny/phynteny.gbk -f {params.pkl}/"$data"_phynteny/phynteny.functions
            done
        else
            touch {output.pkl_func}
        fi
        """

rule annotate_summary_paired:
    input:
        pharokka_func=os.path.join(dir_annotate, "pharokka-pr", "{sample}_1_pharokka", "{sample}_1_pharokka.functions"),
        phold_func=os.path.join(dir_annotate, "phold-pr", "{sample}_1_phold", "{sample}_1_phold.functions"),
        pkl_func=os.path.join(dir_annotate, "phynteny-pr", "{sample}_1_phynteny", "phynteny.functions"),
    params:
        pharokka_func=os.path.join(dir_annotate, "pharokka-pr"),
        phold_func=os.path.join(dir_annotate, "phold-pr"),
        pkl_func=os.path.join(dir_annotate, "phynteny-pr"),
        output=os.path.join(dir_final, "{sample}-pr"),
        ids="{sample}",
    output:
        summary_gbk=os.path.join(dir_final, "{sample}-pr", "{sample}_1_summary.functions")
    localrule: True
    script:
        os.path.join(dir_script, "summary_functions.py")

rule summarize_paired:
    input:
        #this is literally here to make sure there is a file for each sample
        r=os.path.join(dir_fastp, "{sample}_fastp.txt"),
        assembly=os.path.join(dir_megahit, "{sample}-pr", "log"),
        table= os.path.join(dir_genome, "{sample}-pr", "{sample}-genome-candidates.csv"),
        genome=os.path.join(dir_genome, "{sample}-pr", "{sample}.fasta"),
        gbk=os.path.join(dir_annotate, "phynteny-pr", "{sample}_1_phynteny", "phynteny.gbk"),
        plot=os.path.join(dir_annotate, "phynteny-pr", "{sample}_1_phynteny", "plots", "{sample}_1.png"),
        ph_taxa=os.path.join(dir_annotate, "pharokka-pr", "{sample}_1_pharokka", "{sample}_1_top_hits_mash_inphared.tsv"),
        cdden=os.path.join(dir_annotate, "pharokka-pr", "{sample}_1_pharokka", "{sample}_1_length_gc_cds_density.tsv"),
        cds=os.path.join(dir_annotate, "pharokka-pr", "{sample}_1_pharokka", "{sample}_1_cds_functions.tsv"),
        amr =os.path.join(dir_annotate, "pharokka-pr", "{sample}_1_pharokka", "top_hits_card.tsv"),
        vfdb=os.path.join(dir_annotate, "pharokka-pr", "{sample}_1_pharokka", "top_hits_vfdb.tsv"),
        spacers=os.path.join(dir_annotate, "pharokka-pr", "{sample}_1_pharokka",  "{sample}_1_minced_spacers.txt"),
        acr=os.path.join(dir_annotate, "phold-pr", "{sample}_1_phold", "sub_db_tophits", "acr_cds_predictions.tsv"),
        card=os.path.join(dir_annotate, "phold-pr", "{sample}_1_phold", "sub_db_tophits", "card_cds_predictions.tsv"),
        defense=os.path.join(dir_annotate, "phold-pr", "{sample}_1_phold", "sub_db_tophits", "defensefinder_cds_predictions.tsv"),
        vfdb_phold=os.path.join(dir_annotate, "phold-pr", "{sample}_1_phold", "sub_db_tophits", "vfdb_cds_predictions.tsv"),
    output:
        summary=os.path.join(dir_final, "{sample}-pr", "{sample}_summary.txt")
    params:
        #input files for a sample 
        annot=os.path.join(dir_annotate),
        #output files
        genomes= os.path.join(dir_final, "{sample}-pr", "{sample}_genome.fasta"),
        gbks=os.path.join(dir_final, "{sample}-pr", "{sample}.gbk"),
        plots=os.path.join(dir_final, "{sample}-pr", "{sample}_phynteny_plot.png"),
        outdir=os.path.join(dir_final,"{sample}-pr"),
        ID="{sample}",
        seq= "pr"
    localrule: True
    script:
        os.path.join(dir_script, 'summary.py')

rule summarize_annotations_longreads:
    input:
        pharokka=os.path.join(dir_annotate, "pharokka-sr", "{sample}_1_pharokka", "{sample}_1.gbk"),
        phold=os.path.join(dir_annotate, "phold-sr", "{sample}_1_phold", "{sample}_1.gbk"),
        pkl=os.path.join(dir_annotate, "phynteny-sr", "{sample}_1_phynteny", "phynteny.gbk")
    params:
        inputdir=os.path.join(dir_annotate, "{sample}-sr-genomes"),
        pharokka=os.path.join(dir_annotate, "pharokka-sr"),
        phold=os.path.join(dir_annotate, "phold-sr"),
        pkl=os.path.join(dir_annotate, "phynteny-sr")
    output:
        pharokka_func=os.path.join(dir_annotate, "pharokka-sr", "{sample}_1_pharokka", "{sample}_1_pharokka.functions"),
        phold_func=os.path.join(dir_annotate, "phold-sr", "{sample}_1_phold", "{sample}_1_phold.functions"),
        pkl_func=os.path.join(dir_annotate, "phynteny-sr", "{sample}_1_phynteny", "phynteny.functions"),
    resources:
        mem =config['resources']['smalljob']['mem'],
        time = config['resources']['smalljob']['time']
    conda:
        os.path.join(dir_env, "pharokka.yaml")
    shell:
                """
        if [[ -s {input.pharokka} ]] ; then
            for f in {params.inputdir}/*; do 
                data="$(basename "$f" .fasta)"
                genbank_to -g {params.pharokka}/"$data"_pharokka/"$data".gbk -f {params.pharokka}/"$data"_pharokka/"$data"_pharokka.functions
            done
        else
            touch {output.pharokka_func}
        fi

        if [[ -s {input.phold} ]] ; then
            for f in {params.inputdir}/*; do 
                data="$(basename "$f" .fasta)"
                genbank_to -g {params.phold}/"$data"_phold/"$data".gbk -f {params.phold}/"$data"_phold/"$data"_phold.functions
            done
        else
            touch {output.phold_func}
        fi

        if [[ -s {input.pkl} ]] ; then
            for f in {params.inputdir}/*; do 
                data="$(basename "$f" .fasta)"
                genbank_to -g {params.pkl}/"$data"_phynteny/phynteny.gbk -f {params.pkl}/"$data"_phynteny/phynteny.functions
            done
        else
            touch {output.pkl_func}
        fi
        """

rule annotate_summary_longreads:
    input:
        pharokka_func=os.path.join(dir_annotate, "pharokka-sr", "{sample}_1_pharokka", "{sample}_1_pharokka.functions"),
        phold_func=os.path.join(dir_annotate, "phold-sr", "{sample}_1_phold", "{sample}_1_phold.functions"),
        pkl_func=os.path.join(dir_annotate, "phynteny-sr", "{sample}_1_phynteny", "phynteny.functions"),
    params:
        pharokka_func=os.path.join(dir_annotate, "pharokka-sr"),
        phold_func=os.path.join(dir_annotate, "phold-sr"),
        pkl_func=os.path.join(dir_annotate, "phynteny-sr"),
        output=os.path.join(dir_final, "{sample}-sr"),
        ids="{sample}",
    output:
        summary_gbk=os.path.join(dir_final, "{sample}-sr", "{sample}_1_summary.functions")
    localrule: True
    script:
        os.path.join(dir_script, "summary_functions.py")

rule summarize_longread:
    input:
        r=os.path.join(dir_nanopore, "{sample}_filt.txt"),
        assembly=os.path.join(dir_flye, "{sample}-sr", "flye.log"),
        table= os.path.join(dir_genome, "{sample}-sr", "{sample}-genome-candidates.csv"),
        genome=os.path.join(dir_genome, "{sample}-sr", "{sample}.fasta"),
        gbk=os.path.join(dir_annotate, "phynteny-sr", "{sample}_1_phynteny", "phynteny.gbk"),
        plot=os.path.join(dir_annotate, "phynteny-sr", "{sample}_1_phynteny", "plots", "{sample}_1.png"),
        ph_taxa =os.path.join(dir_annotate, "pharokka-sr","{sample}_1_pharokka", "{sample}_1_top_hits_mash_inphared.tsv"),
        cdden=os.path.join(dir_annotate, "pharokka-sr", "{sample}_1_pharokka", "{sample}_1_length_gc_cds_density.tsv"),
        cds=os.path.join(dir_annotate, "pharokka-sr", "{sample}_1_pharokka", "{sample}_1_cds_functions.tsv"),
        amr =os.path.join(dir_annotate, "pharokka-sr", "{sample}_1_pharokka", "top_hits_card.tsv"),
        vfdb=os.path.join(dir_annotate, "pharokka-sr", "{sample}_1_pharokka", "top_hits_vfdb.tsv"),
        spacers=os.path.join(dir_annotate, "pharokka-sr", "{sample}_1_pharokka", "{sample}_1_minced_spacers.txt"),
        acr=os.path.join(dir_annotate, "phold-sr", "{sample}_1_phold","sub_db_tophits", "acr_cds_predictions.tsv"),
        card=os.path.join(dir_annotate,"phold-sr", "{sample}_1_phold","sub_db_tophits", "card_cds_predictions.tsv"),
        defense=os.path.join(dir_annotate, "phold-sr", "{sample}_1_phold","sub_db_tophits", "defensefinder_cds_predictions.tsv"),
        vfdb_phold=os.path.join(dir_annotate, "phold-sr","{sample}_1_phold","sub_db_tophits", "vfdb_cds_predictions.tsv"),
    output:
        summary=os.path.join(dir_final, "{sample}-sr", "{sample}_summary.txt"),
    params:
        #input files for a sample 
        annot=os.path.join(dir_annotate),
        #output files
        genomes= os.path.join(dir_final, "{sample}-sr", "{sample}_genome.fasta"),
        gbks=os.path.join(dir_final, "{sample}-sr", "{sample}.gbk"),
        ID="{sample}",
        plots=os.path.join(dir_final, "{sample}-sr", "{sample}_phynteny_plot.png"),
        outdir=os.path.join(dir_final, "{sample}-sr"),
        seq= "sr"
    localrule: True
    script:
        os.path.join(dir_script, 'summary.py')

