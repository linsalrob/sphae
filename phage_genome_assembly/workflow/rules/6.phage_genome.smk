"""
The resulting assemblies can have multiple contigs, so selecting the phage contigs
"""

rule genomes_megahit:
    input:
        csv = os.path.join(dir.assembly, "{sample}-assembly-stats_megahit.csv")
    output:
        out = os.path.join(dir.genome, "{sample}-pr", "{sample}-genome-candidates.csv")
    conda:
        os.path.join(dir.env, "graph.yaml")
    localrule: True
    log:
        os.path.join(dir.log, "picking_contigs.{sample}.log")
    script:
        os.path.join(dir.script, 'pick_phage_contigs.py')

rule genomes_extract_megahit:
    input:
        contigs = os.path.join(dir.megahit, "{sample}-pr", "final.contigs.fa"),
        csv = os.path.join(dir.genome, "{sample}-pr", "{sample}-genome-candidates.csv"),
    output:
        fasta=os.path.join(dir.genome, "{sample}-pr", "{sample}.fasta")  
    params:
        outdir = os.path.join(dir.genome, "{sample}-pr"),
    conda:
        os.path.join(dir.env, "samtools.yaml")
        
    threads:
        config.resources.smalljob.cpu
    resources:
        mem_mb = config.resources.smalljob.mem,
        time = config.resources.smalljob.time
    log:
        os.path.join(dir.log, "samtools.{sample}.log")
    shell:
        """
        #get the contig or contigs name from the csv file, and run samtools 
        awk -F, 'NR>1 {{print $3}}' {input.csv} > {params.outdir}/phage-genome-contig

        touch {output.fasta}

        #extracting the contigs from the assembly
        for f in `cat {params.outdir}/phage-genome-contig`; do samtools faidx {input.contigs} "$f" >> {output.fasta} ; done 
        """    

# Define a checkpoint rule to check the existence of the output file
checkpoint check_genomes_extract_output:
    input:
        fasta=os.path.join(dir.genome, "{sample}-pr", "{sample}.fasta")        
    output:
        done=os.path.join(dir.genome, "{sample}-pr", "{sample}_genomes_extract_done.txt")
    run:
        import subprocess
        with open(input[0], 'r') as infile:
            if not any(infile):
                # Handle the case where the output file does not exist
                # Running trimnami.smk again but with params set to subsample
                shell("snakemake --use-conda --restart trimnami --config subsample='--subsample' ")    
                
            open(output.done, 'w').close()

rule genomes_flye:
    input:
        csv = os.path.join(dir.assembly, "{sample}-assembly-stats_flye.csv")
    output:
        out =os.path.join(dir.genome, "{sample}-sr", "{sample}-genome-candidates.csv")
    conda:
        os.path.join(dir.env, "graph.yaml")
    localrule: True
    log:
        os.path.join(dir.log, "picking_contigs.{sample}.log")
    script:
        os.path.join(dir.script, 'pick_phage_contigs.py')

rule check_output_length_flye:
    input:
        os.path.join(dir.genome, "{sample}-sr", "{sample}-genome-candidates.csv")
    output:
        os.path.join(dir.genome,"{sample}-sr", "{sample}_length_check.txt")
    localrule: True
    shell:
        """
        output_file="{input}"
        length=$(wc -l <"$output_file")
        echo $length

        if [ $length -ne 2 ]; then
            echo "Entering here"
            #snakemake --snakefile opt.subsampling_long.snakefile --cores all
        elif [ $length -eq 2 ]; then
            echo $length > {output}
        fi
        """

rule genomes_extract_flye:
    input:
        contigs = os.path.join(dir.flye, "{sample}-sr", "assembly.fasta"),
        csv = os.path.join(dir.genome, "{sample}-sr", "{sample}-genome-candidates.csv"),
        lens=os.path.join(dir.genome,"{sample}-sr", "{sample}_length_check.txt")
    output:
        os.path.join(dir.genome, "{sample}-sr", "{sample}.fasta")
    params:
        outdir = os.path.join(dir.genome, "{sample}-sr"),
        sample = "{sample}"
    conda:
        os.path.join(dir.env, "samtools.yaml")
    threads:
        config.resources.smalljob.cpu
    resources:
        mem_mb = config.resources.smalljob.mem,
        time = config.resources.smalljob.time
    log:
        os.path.join(dir.log, "samtools.{sample}.log")
    shell:
        """
        #get the contig or contigs name from the csv file, and run samtools 
        awk -F, 'NR>1 {{print $3}}' {input.csv} > {params.outdir}/phage-genome-contig

        touch {output}
        
        #extracting the contigs from the assembly
        for f in `cat {params.outdir}/phage-genome-contig`; do samtools faidx {input.contigs} "$f" >> {output} ; done 
        """ 

