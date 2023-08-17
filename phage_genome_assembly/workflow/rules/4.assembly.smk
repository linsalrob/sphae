rule flye:
    """Assemble longreads with Flye"""
    input:
        os.path.join(dir.nanopore, "{sample}_single.fastq.gz")
    threads:
        config.resources.bigjob.cpu
    resources:
        mem_mb=config.resources.bigjob.mem,
        time=config.resources.bigjob.time
    output:
        fasta = os.path.join(dir.flye, "{sample}-sr", "assembly.fasta"),
        gfa = os.path.join(dir.flye, "{sample}-sr", "assembly_graph.gfa"),
        path= os.path.join(dir.flye, "{sample}-sr", "assembly_info.txt")
    params:
        out= os.path.join(dir.flye, "{sample}-sr"),
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
        fasta = os.path.join(dir.flye, "{sample}-sr", "assembly.fasta"),
        fastq = os.path.join(dir.nanopore, "{sample}_single.fastq.gz")
    output:
        fasta = os.path.join(dir.flye,"{sample}-sr", "consensus.fasta")
    conda:
        os.path.join(dir.env, "medaka.yaml")
    params:
        model = config.params.medaka,
        dir= directory(os.path.join(dir.flye,"{sample}-sr"))
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

rule spades:
    """Assemble short reads with SPAdes"""
    input:
        r1 = os.path.join(dir.prinseq, "{sample}_R1.fastq.gz"),
        r2 = os.path.join(dir.prinseq, "{sample}_R2.fastq.gz"),
        s =  os.path.join(dir.prinseq, "{sample}_S.fastq.gz"),
    output:
        os.path.join(dir.spades, "{sample}-pr", "contigs.fasta"),
        os.path.join(dir.spades, "{sample}-pr", "contigs.paths"),
        os.path.join(dir.spades, "{sample}-pr", "assembly_graph_with_scaffolds.gfa"),
    params:
        os.path.join(dir.spades, "{sample}-pr")
    log:
        os.path.join(dir.log, "spades.{sample}.log")
    benchmark:
        os.path.join(dir.bench, "spades.{sample}.txt")
    threads:
        config.resources.bigjob.cpu
    resources:
        mem_mb=config.resources.bigjob.mem,
        time=config.resources.bigjob.time
    conda:
        os.path.join(dir.env, "spades.yaml")
    shell:
        """
        if [[ -d {params} ]]
        then
            rm -rf {params}
        fi
        spades.py \
            -1 {input.r1} \
            -2 {input.r2} \
            -s {input.s} \
            -o {params} \
            -t {threads} \
            --careful \
            2> {log}
        """