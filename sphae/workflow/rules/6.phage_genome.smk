"""
The resulting assemblies can have multiple contigs, so selecting the phage contigs
"""

rule genomes_megahit:
    input:
        csv = os.path.join(dir_assembly, "{sample}-assembly-stats_megahit.csv")
    output:
        out = os.path.join(dir_genome, "{sample}-pr", "{sample}-genome-candidates.csv")
    conda:
        os.path.join(dir_env, "graph.yaml")
    localrule: True
    log:
        os.path.join(dir_log, "picking_contigs.{sample}.log")
    script:
        os.path.join(dir_script, 'pick_phage_contigs.py')

rule genomes_extract_megahit:
    input:
        contigs = os.path.join(dir_megahit, "{sample}-pr", "final.contigs.fa"),
        csv = os.path.join(dir_genome, "{sample}-pr", "{sample}-genome-candidates.csv"),
    output:
        fasta=os.path.join(dir_genome, "{sample}-pr", "{sample}.fasta")  
    params:
        outdir = os.path.join(dir_genome, "{sample}-pr"),
    conda:
        os.path.join(dir_env, "samtools.yaml")
    threads:
        config['resources']['smalljob']['cpu']
    resources:
        mem_mb = config['resources']['smalljob']['mem'],
        time = config['resources']['smalljob']['time']
    log:
        os.path.join(dir_log, "samtools.{sample}.log")
    shell:
        """
        #get the contig or contigs name from the csv file, and run samtools 
        awk -F, 'NR>1 {{print $2}}' {input.csv} > {params.outdir}/phage-genome-contig

        touch {output.fasta}

        #extracting the contigs from the assembly
        for f in `cat {params.outdir}/phage-genome-contig`; do samtools faidx {input.contigs} "$f" >> {output.fasta} ; done 
        """    

rule genomes_flye:
    input:
        csv = os.path.join(dir_assembly, "{sample}-assembly-stats_flye.csv")
    output:
        out =os.path.join(dir_genome, "{sample}-sr", "{sample}-genome-candidates.csv")
    conda:
        os.path.join(dir_env, "graph.yaml")
    localrule: True
    log:
        os.path.join(dir_log, "picking_contigs.{sample}.log")
    script:
        os.path.join(dir_script, 'pick_phage_contigs.py')


rule genomes_extract_flye:
    input:
        contigs = os.path.join(dir_flye, "{sample}-sr", "assembly.fasta"),
        csv = os.path.join(dir_genome, "{sample}-sr", "{sample}-genome-candidates.csv"),
    output:
        os.path.join(dir_genome, "{sample}-sr", "{sample}.fasta")
    params:
        outdir = os.path.join(dir_genome, "{sample}-sr"),
        sample = "{sample}"
    conda:
        os.path.join(dir_env, "samtools.yaml")
    threads:
        config['resources']['smalljob']['cpu']
    resources:
        mem_mb = config['resources']['smalljob']['mem'],
        time = config['resources']['smalljob']['time']
    log:
        os.path.join(dir_log, "samtools.{sample}.log")
    shell:
        """
        #get the contig or contigs name from the csv file, and run samtools 
        awk -F, 'NR>1 {{print $2}}' {input.csv} > {params.outdir}/phage-genome-contig

        touch {output}
        
        #extracting the contigs from the assembly
        for f in `cat {params.outdir}/phage-genome-contig`; do samtools faidx {input.contigs} "$f" >> {output} ; done 
        """ 

