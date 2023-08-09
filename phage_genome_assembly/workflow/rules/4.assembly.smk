rule flye:
    """Assemble longreads with Flye"""
    input:
        os.path.join(dir.nanopore, "{sample}.single.fastq.gz")
    threads:
        config.resources.bigjob.cpu
    resources:
        mem_mb=config.resources.bigjob.mem,
        time=config.resources.bigjob.time
    output:
        fasta = os.path.join(dir.flye, "{sample}", "assembly.fasta"),
        gfa = os.path.join(dir.flye, "{sample}", "assembly_graph.gfa"),
        path= os.path.join(dir.flye, "{sample}", "assembly_info.txt")
    params:
        out= os.path.join(dir.flye, "{sample}"),
        model = config.params.flye
    log:
        os.path.join(dir.log, "flye.{sample}.log")
    benchmark:
        os.path.join(dir.bench,"flye.{sample}.txt")
    conda:
        os.path.join(dir.env, "flye.yaml")
    shell:
        """
        flye \
            {params.model} \
            {input} \
            --threads {threads} \
            --out-dir {params.out} \
            2> {log}
        """


rule medaka:
    """Polish longread assembly with medaka"""
    input:
        fasta = os.path.join(dir.flye, "{sample}", "assembly.fasta"),
        fastq = os.path.join(dir.nanopore, "{sample}.single.fastq.gz")
    output:
        fasta = os.path.join(dir.flye,"{sample}", "consensus.fasta")
    conda:
        os.path.join(dir.env, "medaka.yaml")
    params:
        model = config.params.medaka,
        dir= directory(os.path.join(dir.flye,"{sample}"))
    threads:
        config.resources.bigjob.cpu
    resources:
        mem_mb=config.resources.bigjob.mem,
        time=config.resources.bigjob.time
    log:
        os.path.join(dir.log, "medaka.{sample}.log")
    benchmark:
        os.path.join(dir.bench, "medaka.{sample}.txt")
    shell:
        """
        medaka_consensus \
            -i {input.fastq} \
            -d {input.fasta} \
            -o {params.dir} \
            -m {params.model} \
            -t {threads} \
            2> {log}
        """


rule megahit:
    """Assemble short reads with MEGAHIT"""
    input:
        r1 = os.path.join(dir.prinseq, "{sample}_R1.fastq.gz"),
        r2 = os.path.join(dir.prinseq, "{sample}_R2.fastq.gz"),
        s =  os.path.join(dir.prinseq, "{sample}_S.fastq.gz"),
    output:
        os.path.join(dir.megahit, "{sample}", "final.contigs.fa")
    params:
        os.path.join(dir.megahit, "{sample}")
    log:
        os.path.join(dir.log, "megahit.{sample}.log")
    benchmark:
        os.path.join(dir.bench, "megahit.{sample}.txt")
    threads:
        config.resources.bigjob.cpu
    resources:
        mem_mb=config.resources.bigjob.mem,
        time=config.resources.bigjob.time
    conda:
        os.path.join(dir.env, "megahit.yaml")
    shell:
        """
        if [[ -d {params} ]]
        then
            rm -rf {params}
        fi
        megahit \
            -1 {input.r1} \
            -2 {input.r2} \
            -r {input.s} \
            -o {params} \
            -t {threads} \
            2> {log}
        """


rule fastg:
    """Save the MEGAHIT graph"""
    input:
        os.path.join(dir.megahit, "{sample}", "final.contigs.fa")
    output:
        os.path.join(dir.megahit, "{sample}", "final.fastg")
    conda:
        os.path.join(dir.env, "megahit.yaml")
    log:
        os.path.join(dir.log, "fastg.{sample}.log")
    benchmark:
        os.path.join(dir.bench, "fastg.{sample}.txt")
    shell:
        """
        if [[ -s {input} ]]
        then
            kmer=$(head -1 {input} | sed 's/>//' | sed 's/_.*//')
            megahit_toolkit contig2fastg $kmer {input} > {output} 2> {log}
        fi
        """
