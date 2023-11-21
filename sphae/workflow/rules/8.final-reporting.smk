"""
Summarizing the results to one directory
"""

rule summarize_paired:
    input:
        assembly=os.path.join(dir.megahit, "{sample}-pr", "log"),
    output:
        summary=os.path.join(dir.final, "{sample}-pr", "{sample}_summary.txt")
    params:
        genomes= os.path.join(dir.final, "{sample}-pr", "{sample}_genome.fasta"),
        gbks=os.path.join(dir.final, "{sample}-pr", "{sample}.gbk"),
        plots=os.path.join(dir.final, "{sample}-pr", "{sample}_pharokka_plot.png"),
        sample="{sample}",
        reorient=os.path.join(dir.pharokka, "{sample}-pr", "{sample}_dnaapler_reoriented.fasta"),
        genome=os.path.join(dir.genome, "{sample}-pr", "{sample}.fasta"),
        gbk=os.path.join(dir.pharokka, "{sample}-pr", "phynteny", "phynteny.gbk"),
        table= os.path.join(dir.genome, "{sample}-pr", "{sample}-genome-candidates.csv"),
        amr =os.path.join(dir.pharokka, "{sample}-pr", "top_hits_card.tsv"),
        vfdb=os.path.join(dir.pharokka, "{sample}-pr", "top_hits_vfdb.tsv"),
        spacers=os.path.join(dir.pharokka, "{sample}-pr", "{sample}_minced_spacers.txt"),
        plot=os.path.join(dir.pharokka, "{sample}-pr", "{sample}_pharokka_plot.png"),
        ph_taxa =os.path.join(dir.pharokka, "{sample}-pr", "{sample}_top_hits_mash_inphared.tsv"),
        cdden=os.path.join(dir.pharokka, "{sample}-pr", "{sample}_length_gc_cds_density.tsv"),
    localrule: True
    log: 
        os.path.join(dir.log, "final_summary.{sample}.log")
    shell:
        """
        #first checking if the assembly worked
        if (tail -n 1 {input.assembly} | grep 'ALL DONE'); then

            if [ -f "{params.table}" ]; then
                line_count=$(wc -l < "{params.table}")
                if [ "$line_count" -gt 1 ]; then
                    #copying  the final files over
                    if [[ -s {params.reorient} ]]; then
                        cp -r {params.reorient} {params.genomes}
                    else
                        cp -r {params.genome} {params.genomes}
                    fi
                

                    cp -r {params.gbk} {params.gbks}
                    cp -r {params.plot} {params.plots}

                    echo "Sample: $(sed -n '2p' {params.table} | cut -d ',' -f 2)" > {output.summary}
                    
                    line_count=$(wc -l < {params.table})

                    if [ "$line_count" -gt 2 ]; then
                        echo "Genome is incomplete or contaminated, includes mutliple contigs" >> {output.summary}
                    else
                        echo "Length: $(tail -n +2 {params.table} | cut -d ',' -f 13)" >> {output.summary}
                        echo "Coding density: $(tail -n +2 {params.cdden} | cut -f 5 )" >> {output.summary}
                        echo "Circular: $(tail -n +2 {params.table} | cut -d ',' -f 14)" >> {output.summary}
                        echo "Completeness: $(tail -n +2 {params.table} | cut -d ',' -f 21)" >> {output.summary}
                        echo "Contamination: $(tail -n +2 {params.table} | cut -d ',' -f 23)" >> {output.summary}
                        echo "Accession: $(tail -n +2 {params.ph_taxa} | cut -f 2 )" >> {output.summary}
                        echo "Taxa mash hits: $(tail -n +2 {params.ph_taxa} | cut -f 5 )" >> {output.summary}
                        echo "Taxa_Name: $(tail -n +2 {params.ph_taxa} | cut -f 6 )" >> {output.summary}
                    fi

                    if grep -q "int" {params.gbk}; then
                        echo "Integrase found" >> {output.summary}
                    else
                        echo "No integrase" >> {output.summary}
                    fi

                    if [[ $(wc -l < "{params.amr}") -eq 1 ]]; then
                        echo "No AMR genes" >> {output.summary}
                    else
                        echo "AMR genes found" >> {output.summary}
                    fi

                    if [[ $(wc -l < "{params.vfdb}") -eq 1 ]]; then
                        echo "No virulence factor genes" >> {output.summary}
                    else
                        echo "Virulence factor genes found" >> {output.summary}
                    fi

                    if [[ $(wc -l < "{params.spacers}") -gt 2 ]]; then
                        echo "CRISPR spacers found" >> {output.summary}
                    else
                        echo "No CRISPR spacers found" >> {output.summary}
                    fi
                else
                    #this is if the assembly is too fragmented
                    echo "Sample: {params.sample}" > {output.summary}
                    echo "Assembly likely too fragmented or none of the contigs were assigned viral" >> {output.summary}
                fi
            
            else
                #this is if the assembly is too fragmented
                echo "Sample: {params.sample}" > {output.summary}
                echo "Assembly likely too fragmented or none of the contigs were assigned viral" >> {output.summary}
            fi
        else
            #this is if the assembly failed then
            echo "Sample: {params.sample}" > {output.summary}
            echo "Failed during assembly" >> {output.summary}
        fi
        """

rule summarize_longread:
    input:
        assembly=os.path.join(dir.flye, "{sample}-sr", "flye.log"),
    output:
        summary=os.path.join(dir.final, "{sample}-sr", "{sample}_summary.txt")
    params:
        genomes= os.path.join(dir.final, "{sample}-sr", "{sample}_genome.fasta"),
        gbks=os.path.join(dir.final, "{sample}-sr", "{sample}.gbk"),
        plots=os.path.join(dir.final, "{sample}-sr", "{sample}_pharokka_plot.png"),
        sample="{sample}",
        reorient=os.path.join(dir.pharokka, "{sample}-sr", "{sample}_dnaapler_reoriented.fasta"),
        genome=os.path.join(dir.genome, "{sample}-sr", "{sample}.fasta"),
        gbk=os.path.join(dir.pharokka, "{sample}-sr", "phynteny", "phynteny.gbk"),
        table= os.path.join(dir.genome, "{sample}-sr", "{sample}-genome-candidates.csv"),
        amr =os.path.join(dir.pharokka, "{sample}-sr", "top_hits_card.tsv"),
        vfdb=os.path.join(dir.pharokka, "{sample}-sr", "top_hits_vfdb.tsv"),
        plot=os.path.join(dir.pharokka, "{sample}-sr", "{sample}_pharokka_plot.png"),
        spacers=os.path.join(dir.pharokka, "{sample}-sr", "{sample}_minced_spacers.txt"),
        ph_taxa =os.path.join(dir.pharokka, "{sample}-sr", "{sample}_top_hits_mash_inphared.tsv"),
        cdden=os.path.join(dir.pharokka, "{sample}-sr", "{sample}_length_gc_cds_density.tsv"),
    localrule: True
    log: 
        os.path.join(dir.log, "final_summary.{sample}.log")
    shell:
        """
        #first checking if the assembly worked
        if tail -n 1 {input.assembly} | grep 'INFO: Final assembly:' ; then

            if [[ -f "{params.table}" ]]; then
                line_count=$(wc -l < "{params.table}")
                if [ "$line_count" -gt 1 ]; then
                    #copying  the final files over
                    if [[ -s {params.reorient} ]]; then
                        cp -r {params.reorient} {params.genomes}
                    else
                        cp -r {params.genome} {params.genomes}
                    fi
                

                    cp -r {params.gbk} {params.gbks}
                    cp -r {params.plot} {params.plots}

                    echo "Sample: $(sed -n '2p' {params.table} | cut -d ',' -f 2)" > {output.summary}
                    
                    line_count=$(wc -l < {params.table})

                    if [ "$line_count" -gt 2 ]; then
                        echo "Genome is incomplete or contaminated, includes mutliple contigs" >> {output.summary}
                    else
                        echo "Length: $(tail -n +2 {params.table} | cut -d ',' -f 13)" >> {output.summary}
                        echo "Coding density: $(tail -n +2 {params.cdden} | cut -f 5 )" >> {output.summary}
                        echo "Circular: $(tail -n +2 {params.table} | cut -d ',' -f 14)" >> {output.summary}
                        echo "Completeness: $(tail -n +2 {params.table} | cut -d ',' -f 21)" >> {output.summary}
                        echo "Contamination: $(tail -n +2 {params.table} | cut -d ',' -f 23)" >> {output.summary}
                        echo "Accession: $(tail -n +2 {params.ph_taxa} | cut -f 2 )" >> {output.summary}
                        echo "Taxa mash hits: $(tail -n +2 {params.ph_taxa} | cut -f 5 )" >> {output.summary}
                        echo "Taxa_Name: $(tail -n +2 {params.ph_taxa} | cut -f 6 )" >> {output.summary}
                    fi

                    if grep -q "int" {params.gbk}; then
                        echo "Integrase found" >> {output.summary}
                    else
                        echo "No integrase" >> {output.summary}
                    fi

                    if [[ $(wc -l < "{params.amr}") -eq 1 ]]; then
                        echo "No AMR genes" >> {output.summary}
                    else
                        echo "AMR genes found" >> {output.summary}
                    fi

                    if [[ $(wc -l < "{params.vfdb}") -eq 1 ]]; then
                        echo "No virulence factor genes" >> {output.summary}
                    else
                        echo "Virulence factor genes found" >> {output.summary}
                    fi

                    if [[ $(wc -l < "{params.spacers}") -gt 2 ]]; then
                        echo "CRISPR spacers found" >> {output.summary}
                    else
                        echo "No CRISPR spacers found" >> {output.summary}
                    fi
                else
                    #this is if the assembly is too fragmented
                    echo "Sample: {params.sample}" > {output.summary}
                    echo "Assembly likely too fragmented or none of the contigs were assigned viral" >> {output.summary}
                fi
            
            else
                #this is if the assembly is too fragmented
                echo "Sample: {params.sample}" > {output.summary}
                echo "Assembly likely too fragmented or none of the contigs were assigned viral" >> {output.summary}
            fi
        
        else
            #this is if the assembly failed then
            echo "Sample: {params.sample}" > {output.summary}
            echo "Failed during assembly" >> {output.summary}
        fi
        """
