"""
Assembly rules
Illumina paired end reads - Megahit
Nanopore reads - Flye
"""

rule flye:
    """Assemble longreads with Flye"""
    input:
        os.path.join(dir_nanopore, "{sample}_filt.fastq.gz"),
    output:
        fasta = os.path.join(dir_flye, "{sample}-sr", "assembly.fasta"),
        gfa = os.path.join(dir_flye, "{sample}-sr", "assembly_graph.gfa"),
        path= os.path.join(dir_flye, "{sample}-sr", "assembly_info.txt"),
        log=os.path.join(dir_flye, "{sample}-sr", "flye.log")
    params:
        out= os.path.join(dir_flye, "{sample}-sr"),
        model = config['params']['flye'],
        g = config['params']['genomeSize']
    log:
        os.path.join(dir_log, "flye.{sample}.log")
    conda:
        os.path.join(dir_env, "flye.yaml")
    threads:
        config['resources']['bigjob']['cpu']
    resources:
        mem_mb=config['resources']['bigjob']['mem'],
        time=config['resources']['bigjob']['time']
    shell:
        """
        if flye \
            {params.model} \
            {input} \
            --threads {threads}  \
            --asm-coverage 50 \
            --genome-size {params.g} \
            --out-dir {params.out} \
            2> {log}; then
                touch {output.fasta}
                touch {output.gfa}
                touch {output.path}
                touch {output.log}
            else
                touch {output.fasta}
                touch {output.gfa}
                touch {output.path}
                touch {output.log}
        fi
        """

rule medaka:
    """Polish longread assembly with medaka"""
    input:
        fasta = os.path.join(dir_flye, "{sample}-sr", "assembly.fasta"),
        fastq = os.path.join(dir_nanopore, "{sample}_filt.fastq.gz"),
    output:
        fasta = os.path.join(dir_flye,"{sample}-sr", "consensus.fasta")
    conda:
        os.path.join(dir_env, "medaka.yaml")
    params:
        model = config['params']['medaka'],
        dir= directory(os.path.join(dir_flye,"{sample}-sr"))
    threads:
        config['resources']['bigjob']['cpu']
    resources:
        mem_mb=config['resources']['bigjob']['mem'],
        time=config['resources']['bigjob']['time']
    log:
        os.path.join(dir_log, "medaka.{sample}.log")
    shell:
        """
        if [[ -s {input.fasta} ]] ; then
            medaka_consensus \
                -i {input.fastq} \
                -d {input.fasta} \
                -o {params.dir} \
                -m {params.model} \
                -t {threads} \
                2> {log}
            touch {output.fasta}
        else
            touch {output.fasta}
        fi
        """


rule megahit:
    """Assemble short reads with MEGAHIT"""
    input:
        r1 = os.path.join(dir_fastp, "{sample}_subsampled_R1.fastq.gz"),
        r2 = os.path.join(dir_fastp, "{sample}_subsampled_R2.fastq.gz")
    output:
        contigs=os.path.join(dir_megahit, "{sample}-pr", "final.contigs.fa"),
        log=os.path.join(dir_megahit, "{sample}-pr", "log")
    params:
        os.path.join(dir_megahit, "{sample}-pr")
    log:
        os.path.join(dir_log, "megahit.{sample}.log")
    threads:
        config['resources']['bigjob']['cpu']
    resources:
        mem_mb=config['resources']['bigjob']['mem'],
        time=config['resources']['bigjob']['time']
    conda:
        os.path.join(dir_env, "megahit.yaml")
    shell:
        """
        if megahit \
            -1 {input.r1} \
            -2 {input.r2} \
            -o {params} \
            -t {threads} -f \
            2> {log}; then
            touch {output.contigs}
            touch {output.log}
        else
            touch {output.contigs}
            touch {output.log}
        fi
        """

rule fastg:
    """Save the MEGAHIT graph"""
    input:
        os.path.join(dir_megahit, "{sample}-pr", "final.contigs.fa")
    output:
        fastg=os.path.join(dir_megahit, "{sample}-pr", "final.fastg"),
        graph=os.path.join(dir_megahit, "{sample}-pr", "final.gfa")
    conda:
        os.path.join(dir_env, "megahit.yaml")
    log:
        os.path.join(dir_log, "fastg.{sample}.log")
    shell:
        """
        if [[ -s {input} ]] ; then
            kmer=$(head -1 {input} | sed 's/>//' | sed 's/_.*//')
            megahit_toolkit contig2fastg $kmer {input} > {output.fastg} 2> {log}
            Bandage reduce {output.fastg} {output.graph} 2>> {log}
            touch {output.fastg}
            touch {output.graph}
        else
            touch {output.fastg}
            touch {output.graph}
        fi
        """
