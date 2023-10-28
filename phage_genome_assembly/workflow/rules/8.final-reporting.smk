"""
Summarizing the results to one directory
"""

rule summarize_paired:
    input:
        genomes= os.path.join(dir.genome, "{sample}-pr", "{sample}_genome.fasta"),
        gbk=os.path.join(dir.pharokka, "{sample}-pr", "phynteny", "phynteny.gbk"),
        table= os.path.join(dir.genome, "{sample}-pr", "{sample}-genome-candidates.csv"),
        amr =os.path.join(dir.pharokka, "{sample}-pr", "top_hits_card.tsv"),
        vfdb=os.path.join(dir.pharokka, "{sample}-pr", "top_hits_vfdb.tsv"),
        spacers=os.path.join(dir.pharokka, "{sample}-pr", "{sample}_minced_spacers.txt"),
        plot=os.path.join(dir.pharokka, "{sample}-pr", "{sample}_pharokka_plot.png"),
        ph_taxa =os.path.join(dir.pharokka, "{sample}-pr", "{sample}_top_hits_mash_inphared.tsv"),
        taxa=os.path.join(dir.taxa, "{sample}-pr", "Summary_taxonomy.tsv")
    output:
        summary=os.path.join(dir.final, "{sample}-pr", "{sample}_summary.txt")
    params:
        genomes= os.path.join(dir.final, "{sample}-pr", "{sample}_genome.fasta"),
        gbk=os.path.join(dir.final, "{sample}-pr", "{sample}.gbk"),
        plot=os.path.join(dir.final, "{sample}-pr", "{sample}_pharokka_plot.png"),
    localrule: True
    log: 
        os.path.join(dir.log, "final_summary.{sample}.log")
    shell:
        """
        #copying the final files over
        cp -r {input.genomes} {params.genomes}
        cp -r {input.gbk} {params.gbk}
        cp -r {input.plot} {params.plot}

        echo "Sample: $(tail -n +2 {input.table} | cut -d ',' -f 2)" > {output.summary}
        echo "Length: $(tail -n +2 {input.table} | cut -d ',' -f 13)" >> {output.summary}
        echo "Circular: $(tail -n +2 {input.table} | cut -d ',' -f 14)" >> {output.summary}
        echo "Completeness: $(tail -n +2 {input.table} | cut -d ',' -f 31)" >> {output.summary}
        echo "Contamination: $(tail -n +2 {input.table} | cut -d ',' -f 33)" >> {output.summary}
        echo "Accession: $(tail -n +2 {input.ph_taxa} | cut -f 2 )" >> {output.summary}
        echo "Taxa_Name: $(tail -n +2 {input.ph_taxa} | cut -f 6 )" >> {output.summary}

        if grep -q "int" {input.gbk}; then
            echo "Integrase found" >> {output.summary}
        else
            echo "No integrase" >> {output.summary}
        fi

        if [[ $(wc -l < "{input.amr}") -eq 1 ]]; then
            echo "No AMR genes" >> {output.summary}
        else
            echo "AMR genes found" >> {output.summary}
        fi

        if [[ $(wc -l < "{input.vfdb}") -eq 1 ]]; then
            echo "No virulence factor genes" >> {output.summary}
        else
            echo "Virulence factor genes found" >> {output.summary}
        fi

        if [[ $(wc -l < "{input.spacers}") -eq 1 ]]; then
            echo "No CRISPR spacers found" >> {output.summary}
        else
            echo "CRISPR spacers found" >> {output.summary}
        fi
        """

rule summarize_longread:
    input:
        genomes= os.path.join(dir.genome, "{sample}-sr", "{sample}_genome.fasta"),
        gbk=os.path.join(dir.pharokka, "{sample}-sr", "phynteny", "phynteny.gbk"),
        table= os.path.join(dir.genome, "{sample}-sr", "{sample}-genome-candidates.csv"),
        amr =os.path.join(dir.pharokka, "{sample}-sr", "top_hits_card.tsv"),
        vfdb=os.path.join(dir.pharokka, "{sample}-sr", "top_hits_vfdb.tsv"),
        plot=os.path.join(dir.pharokka, "{sample}-sr", "{sample}_pharokka_plot.png"),
        spacers=os.path.join(dir.pharokka, "{sample}-sr", "{sample}_minced_spacers.txt"),
        ph_taxa =os.path.join(dir.pharokka, "{sample}-sr", "{sample}_top_hits_mash_inphared.tsv"),
        taxa=os.path.join(dir.taxa, "{sample}-sr", "Summary_taxonomy.tsv")
    output:
        summary=os.path.join(dir.final, "{sample}-sr", "{sample}_summary.txt")
    params:
        genomes= os.path.join(dir.final, "{sample}-sr", "{sample}_genome.fasta"),
        gbk=os.path.join(dir.final, "{sample}-sr", "{sample}.gbk"),
        plot=os.path.join(dir.final, "{sample}-sr", "{sample}_pharokka_plot.png"),
    localrule: True
    log: 
        os.path.join(dir.log, "final_summary.{sample}.log")
    shell:
        """
        #copying the final files over
        cp -r {input.genomes} {params.genomes}
        cp -r {input.gbk} {params.gbk}
        cp -r {input.plot} {params.plot}
  
        echo "Sample: $(tail -n +2 {input.table} | cut -d ',' -f 2)" > {output.summary}
        echo "Length: $(tail -n +2 {input.table} | cut -d ',' -f 13)" >> {output.summary}
        echo "Circular: $(tail -n +2 {input.table} | cut -d ',' -f 14)" >> {output.summary}
        echo "Completeness: $(tail -n +2 {input.table} | cut -d ',' -f 31)" >> {output.summary}
        echo "Contamination: $(tail -n +2 {input.table} | cut -d ',' -f 33)" >> {output.summary}
        echo "Accession: $(tail -n +2 {input.ph_taxa} | cut -f 2 )" >> {output.summary}
        echo "Taxa_Name: $(tail -n +2 {input.ph_taxa} | cut -f 6 )" >> {output.summary}

        if grep -q "int" {input.gbk}; then
            echo "Integrase found" >> {output.summary}
        else
            echo "No integrase" >> {output.summary}
        fi

        if [[ $(wc -l < "{input.amr}") -eq 1 ]]; then
            echo "No AMR genes" >> {output.summary}
        else
            echo "AMR genes found" >> {output.summary}
        fi

        if [[ $(wc -l < "{input.vfdb}") -eq 1 ]]; then
            echo "No virulence factor genes" >> {output.summary}
        else
            echo "Virulence factor genes found" >> {output.summary}
        fi

        if [[ $(wc -l < "{input.spacers}") -eq 1 ]]; then
            echo "No CRISPR spacers found" >> {output.summary}
        else
            echo "CRISPR spacers found" >> {output.summary}
        fi
        """