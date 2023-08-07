"""
The resulting assemblies can have multiple contigs, so selecting the phage contigs
"""

rule genomes_megahit:
    input:
        csv = os.path.join(dir.assembly, "{sample}-assembly-stats_megahit.csv")
    output:
        out =os.path.join(dir.genome, "{sample}", "{sample}-contig-candidates.csv")
    conda:
        os.path.join(dir.env, "graph.yaml")
    localrule: True
    log:
        os.path.join(dir.log, "picking_contigs.{sample}.log")
    script:
        os.path.join(dir.script, 'picking_contigs.py')

rule genomes_extract_megahit:
    input:
        contigs = os.path.join(dir.megahit, "{sample}", "final.contigs.fa"),
        csv = os.path.join(dir.genome, "{sample}", "{sample}-contig-candidates.csv")
    output:
        os.path.join(dir.genome, "{sample}", "{sample}.fasta")
    params:
        outdir = os.path.join(dir.genome, "{sample}"),
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
        cut -f 3 -d "," {input.csv} > tmp

        touch {output}

        #extracting the contigs from the assembly
        for f in `cat tmp` do; samtools faidx {input.contigs} '$f' > {params.outdir}/'$f'.fasta ; done 

        #removing the tmp file
        rm -rf tmp
        """    