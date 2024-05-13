"""
Running pharokka for anntoation 
"""
rule rename_contigs_megahit:
    input:
        fin=os.path.join(dir.genome, "{sample}-pr", "{sample}.fasta"),
    params:
        s ="{sample}"
    output:
        out=os.path.join(dir.genome, "{sample}-pr", "{sample}_genome.fasta"),
        csv=os.path.join(dir.genome, "{sample}-pr", "{sample}_temp.csv")
    localrule: True
    log:
        os.path.join(dir.log, "rename-contigs.{sample}.log")
    script:
        os.path.join(dir.script, 'rename_genomes.py')

rule rename_contigs_flye:
    input:
        fin=os.path.join(dir.genome, "{sample}-sr", "{sample}.fasta"),
    params:
        s ="{sample}"
    output:
        out=os.path.join(dir.genome, "{sample}-sr", "{sample}_genome.fasta"),
        csv=os.path.join(dir.genome, "{sample}-sr", "{sample}_temp.csv")
    localrule: True
    log:
        os.path.join(dir.log, "rename-contigs.{sample}.log")
    script:
        os.path.join(dir.script, 'rename_genomes.py')
