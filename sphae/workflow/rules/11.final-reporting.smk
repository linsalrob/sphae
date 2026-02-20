"""
Summarizing the results to one directory
"""
rule summarize_annotations_paired:
    input: 
        pharokka=os.path.join(dir_annot, "pharokka-pr", "{sample}_1_pharokka", "{sample}_1.gbk"),
        phold=os.path.join(dir_annot, "phold-pr", "{sample}_1_phold", "{sample}_1.gbk"),
        pkl=os.path.join(dir_annot, "phynteny-pr", "{sample}_1_phynteny", "phynteny.gbk")
    params:
        inputdir=os.path.join(dir_annot, "{sample}-pr-genomes"),
        pharokka=os.path.join(dir_annot, "pharokka-pr"),
        phold=os.path.join(dir_annot, "phold-pr"),
        pkl=os.path.join(dir_annot, "phynteny-pr"),
    output:
        pharokka_func=os.path.join(dir_annot, "pharokka-pr", "{sample}_1_pharokka", "{sample}_1_pharokka.functions"),
        phold_func=os.path.join(dir_annot, "phold-pr", "{sample}_1_phold", "{sample}_1_phold.functions"),
        pkl_func=os.path.join(dir_annot, "phynteny-pr", "{sample}_1_phynteny", "phynteny.functions"),
    resources:
        mem_mb =config['resources']['smalljob']['mem_mb'],
        runtime = config['resources']['smalljob']['runtime']
    conda:
        os.path.join(dir_env, "phold.yaml")
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
        pharokka_func=os.path.join(dir_annot, "pharokka-pr", "{sample}_1_pharokka", "{sample}_1_pharokka.functions"),
        phold_func=os.path.join(dir_annot, "phold-pr", "{sample}_1_phold", "{sample}_1_phold.functions"),
        pkl_func=os.path.join(dir_annot, "phynteny-pr", "{sample}_1_phynteny", "phynteny.functions"),
    params:
        pharokka_func=os.path.join(dir_annot, "pharokka-pr"),
        phold_func=os.path.join(dir_annot, "phold-pr"),
        pkl_func=os.path.join(dir_annot, "phynteny-pr"),
        tmp=os.path.join(dir_annot, "phynteny-pr", "tmp"),
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
        gbk=os.path.join(dir_annot, "phynteny-pr", "{sample}_1_phynteny", "phynteny.gbk"),
        plot=os.path.join(dir_annot, "phynteny-pr", "{sample}_1_phynteny", "plots", "{sample}_1.png"),
        ph_taxa=os.path.join(dir_annot, "pharokka-pr", "{sample}_1_pharokka", "{sample}_1_top_hits_mash_inphared.tsv"),
        cdden=os.path.join(dir_annot, "pharokka-pr", "{sample}_1_pharokka", "{sample}_1_length_gc_cds_density.tsv"),
        cds=os.path.join(dir_annot, "pharokka-pr", "{sample}_1_pharokka", "{sample}_1_cds_functions.tsv"),
        amr =os.path.join(dir_annot, "pharokka-pr", "{sample}_1_pharokka", "top_hits_card.tsv"),
        vfdb=os.path.join(dir_annot, "pharokka-pr", "{sample}_1_pharokka", "top_hits_vfdb.tsv"),
        spacers=os.path.join(dir_annot, "pharokka-pr", "{sample}_1_pharokka",  "{sample}_1_minced_spacers.txt"),
        acr=os.path.join(dir_annot, "phold-pr", "{sample}_1_phold", "sub_db_tophits", "acr_cds_predictions.tsv"),
        card=os.path.join(dir_annot, "phold-pr", "{sample}_1_phold", "sub_db_tophits", "card_cds_predictions.tsv"),
        defense=os.path.join(dir_annot, "phold-pr", "{sample}_1_phold", "sub_db_tophits", "defensefinder_cds_predictions.tsv"),
        vfdb_phold=os.path.join(dir_annot, "phold-pr", "{sample}_1_phold", "sub_db_tophits", "vfdb_cds_predictions.tsv"),
    output:
        summary=os.path.join(dir_final, "{sample}-pr", "{sample}_summary.txt")
    params:
        #input files for a sample 
        annot=os.path.join(dir_annot),
        #output files
        genomes= os.path.join(dir_final, "{sample}-pr", "{sample}_genome.fasta"),
        gbks=os.path.join(dir_final, "{sample}-pr", "{sample}.gbk"),
        plots=os.path.join(dir_final, "{sample}-pr", "{sample}_phynteny_plot.png"),
        outdir=os.path.join(dir_final,"{sample}-pr"),
        ID="{sample}",
        seq= "pr"
    script:
        os.path.join(dir_script, 'summary.py')

rule copy_accessory_pr:
    input:
        card=os.path.join(dir_annot, "phold-pr", "{sample}_1_phold", "sub_db_tophits", "card_cds_predictions.tsv"),
        defense=os.path.join(dir_annot, "phold-pr", "{sample}_1_phold", "sub_db_tophits", "defensefinder_cds_predictions.tsv"),
        vfdb_phold=os.path.join(dir_annot, "phold-pr", "{sample}_1_phold", "sub_db_tophits", "vfdb_cds_predictions.tsv"),
    output:
        pseudo=os.path.join(dir_final, "{sample}-pr", "{sample}_tmp"),
    params:
        outdir=os.path.join(dir_final, "{sample}-pr"),
        ptv=os.path.join(dir_phageterm, "{sample}_pr_phageterm"),
        indir=os.path.join(dir_annot, "phold-pr"),
        s="{sample}"
    shell:
        """
        # Loop over ALL genomes for this sample: sample_1_phold, sample_2_phold, ...
        for f in {params.indir}/{params.s}_*_phold; do
            # skip if nothing matches
            [ -e "$f" ] || continue

            # e.g. basename "sample_1_phold" -> "sample_1_phold"
            # then strip the "_phold" suffix -> "sample_1"
            base=$(basename "$f" _phold)

            cp "$f/sub_db_tophits/card_cds_predictions.tsv" \
               "{params.outdir}/"$base"_phold_amr.tsv"

            cp "$f/sub_db_tophits/defensefinder_cds_predictions.tsv" \
               "{params.outdir}/"$base"_phold_defense.tsv"

            cp "$f/sub_db_tophits/vfdb_cds_predictions.tsv" \
               "{params.outdir}/"$base"_phold_vfdb.tsv"
        done

        touch {output.pseudo}
        mv {params.ptv} {params.outdir}/.
        """
    
rule summarize_annotations_longreads:
    input:
        pharokka=os.path.join(dir_annot, "pharokka-sr", "{sample}_1_pharokka", "{sample}_1.gbk"),
        phold=os.path.join(dir_annot, "phold-sr", "{sample}_1_phold", "{sample}_1.gbk"),
        pkl=os.path.join(dir_annot, "phynteny-sr", "{sample}_1_phynteny", "phynteny.gbk")
    params:
        inputdir=os.path.join(dir_annot, "{sample}-sr-genomes"),
        pharokka=os.path.join(dir_annot, "pharokka-sr"),
        phold=os.path.join(dir_annot, "phold-sr"),
        pkl=os.path.join(dir_annot, "phynteny-sr")
    output:
        pharokka_func=os.path.join(dir_annot, "pharokka-sr", "{sample}_1_pharokka", "{sample}_1_pharokka.functions"),
        phold_func=os.path.join(dir_annot, "phold-sr", "{sample}_1_phold", "{sample}_1_phold.functions"),
        pkl_func=os.path.join(dir_annot, "phynteny-sr", "{sample}_1_phynteny", "phynteny.functions"),
    resources:
        mem_mb =config['resources']['smalljob']['mem_mb'],
        runtime = config['resources']['smalljob']['runtime']
    conda:
        os.path.join(dir_env, "phold.yaml")
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
        pharokka_func=os.path.join(dir_annot, "pharokka-sr", "{sample}_1_pharokka", "{sample}_1_pharokka.functions"),
        phold_func=os.path.join(dir_annot, "phold-sr", "{sample}_1_phold", "{sample}_1_phold.functions"),
        pkl_func=os.path.join(dir_annot, "phynteny-sr", "{sample}_1_phynteny", "phynteny.functions"),
    params:
        pharokka_func=os.path.join(dir_annot, "pharokka-sr"),
        phold_func=os.path.join(dir_annot, "phold-sr"),
        pkl_func=os.path.join(dir_annot, "phynteny-sr"),
        output=os.path.join(dir_final, "{sample}-sr"),
        tmp=os.path.join(dir_annot, "phynteny-sr", "tmp"),
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
        gbk=os.path.join(dir_annot, "phynteny-sr", "{sample}_1_phynteny", "phynteny.gbk"),
        plot=os.path.join(dir_annot, "phynteny-sr", "{sample}_1_phynteny", "plots", "{sample}_1.png"),
        ph_taxa =os.path.join(dir_annot, "pharokka-sr","{sample}_1_pharokka", "{sample}_1_top_hits_mash_inphared.tsv"),
        cdden=os.path.join(dir_annot, "pharokka-sr", "{sample}_1_pharokka", "{sample}_1_length_gc_cds_density.tsv"),
        cds=os.path.join(dir_annot, "pharokka-sr", "{sample}_1_pharokka", "{sample}_1_cds_functions.tsv"),
        amr =os.path.join(dir_annot, "pharokka-sr", "{sample}_1_pharokka", "top_hits_card.tsv"),
        vfdb=os.path.join(dir_annot, "pharokka-sr", "{sample}_1_pharokka", "top_hits_vfdb.tsv"),
        spacers=os.path.join(dir_annot, "pharokka-sr", "{sample}_1_pharokka", "{sample}_1_minced_spacers.txt"),
        acr=os.path.join(dir_annot, "phold-sr", "{sample}_1_phold","sub_db_tophits", "acr_cds_predictions.tsv"),
        card=os.path.join(dir_annot,"phold-sr", "{sample}_1_phold","sub_db_tophits", "card_cds_predictions.tsv"),
        defense=os.path.join(dir_annot, "phold-sr", "{sample}_1_phold","sub_db_tophits", "defensefinder_cds_predictions.tsv"),
        vfdb_phold=os.path.join(dir_annot, "phold-sr","{sample}_1_phold","sub_db_tophits", "vfdb_cds_predictions.tsv"),
    output:
        summary=os.path.join(dir_final, "{sample}-sr", "{sample}_summary.txt"),
    params:
        #input files for a sample 
        annot=os.path.join(dir_annot),
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

rule copy_accessory_sr:
    input:
        card=os.path.join(dir_annot, "phold-sr", "{sample}_1_phold", "sub_db_tophits", "card_cds_predictions.tsv"),
        defense=os.path.join(dir_annot, "phold-sr", "{sample}_1_phold", "sub_db_tophits", "defensefinder_cds_predictions.tsv"),
        vfdb_phold=os.path.join(dir_annot, "phold-sr", "{sample}_1_phold", "sub_db_tophits", "vfdb_cds_predictions.tsv"),
    output:
        pseudo=os.path.join(dir_final, "{sample}-sr", "{sample}_tmp"),
    params:
        outdir=os.path.join(dir_final, "{sample}-sr"),
        ptv=os.path.join(dir_phageterm, "{sample}_sr_phageterm"),
        indir=os.path.join(dir_annot, "phold-sr"),
        s="{sample}"
    shell:
        """
        # Loop over ALL genomes for this sample: sample_1_phold, sample_2_phold, ...
        rm -rf {params.outdir}/card_cds_predictions.tsv {params.outdir}/defensefinder_cds_predictions.tsv {params.outdir}/vfdb_cds_predictions.tsv
        for f in {params.indir}/{params.s}_*_phold; do
            # skip if nothing matches
            [ -e "$f" ] || continue

            # e.g. basename "sample_1_phold" -> "sample_1_phold"
            # then strip the "_phold" suffix -> "sample_1"
            base=$(basename "$f" _phold)

            cp "$f/sub_db_tophits/card_cds_predictions.tsv" \
               "{params.outdir}/"$base"_phold_amr.tsv"

            cp "$f/sub_db_tophits/defensefinder_cds_predictions.tsv" \
               "{params.outdir}/"$base"_phold_defense.tsv"

            cp "$f/sub_db_tophits/vfdb_cds_predictions.tsv" \
               "{params.outdir}/"$base"_phold_vfdb.tsv"
        done
        mv {params.ptv} {params.outdir}/.
        touch {output.pseudo}
        """

