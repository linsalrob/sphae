from pathlib import Path
import shutil
import pandas as pd

def copy_files(input_files, params):
    shutil.copy(input_files['genome'], params['genomes'])
    shutil.copy(input_files['gbk'], params['gbks'])
    shutil.copy(input_files['plot'], params['plots'])

def generate_summary(input_files, output_summary, params):
    with open(output_summary, 'w') as summary:
        summary.write(f"Sample: {params['sample']}\n")
        if Path(input_files['table']).exists():
            line_count = sum(1 for line in open(input_files['table']))
            if line_count > 1:
                with open(input_files['table'], 'r') as table_file:
                    lines = table_file.readlines()
                    #print (lines[1])
                    summary.write(f"Length: {lines[1].split(',')[2]}\n")
                    summary.write(f"Circular: {lines[1].split(',')[3]}\n")
                    summary.write(f"Graph connections: {lines[1].split(',')[4]}\n")
                    summary.write(f"Completeness: {lines[1].split(',')[20]}\n")
                    summary.write(f"Contamination: {lines[1].split(',')[22]}\n")
                
                with open(input_files['taxa'], 'r') as taxa_file:
                    tax = pd.read_csv(taxa_file, sep='\t')
                    summary.write("Taxa name (Matching hashes):\t")
                    for index, row in tax.iterrows():
                        summary.write(f"{row['Description']}\t{row['mash_matching_hashes']}\n")

                with open(input_files['cds'], 'r') as cds:
                    cds_df=pd.read_csv(cds, sep='\t')
                    cds_data = cds_df[cds_df['Description'] == 'CDS']
                    count_value = cds_data['Count'].values[0]
                    summary.write(f"Number of CDS: {count_value}\n")
                
                with open(input_files['cdden'], 'r') as cdden:
                    cdn=pd.read_csv(cdden, sep='\t')
                    gc_percent = cdn['gc_perc'].values[0]
                    coding_density = cdn['cds_coding_density'].values[0]
                    summary.write(f"GC percent: {gc_percent}%\n")
                    summary.write(f"Coding density: {coding_density}\n")

                if 'integra' in open(input_files['gbk']).read():
                    summary.write("Integrase found, below is the gene name found\n")
                    with open(input_files['gbk'], 'r') as gbk_file:
                        for line in gbk_file:
                            if 'integra' in line:
                                summary.write(line)
                else:
                    summary.write("No integrase\n")
                
                if 'recombinase' in open(input_files['gbk']).read():
                    summary.write("Recombinase found, below is the gene name found\n")
                    with open(input_files['gbk'], 'r') as gbk_file:
                        for line in gbk_file:
                            if 'recombinase' in line:
                                summary.write(line)
                else:
                    summary.write("No Recombinase\n")

                if 'transposase' in open(input_files['gbk']).read():
                    summary.write("Transposase found, below is the gene name found\n")
                    with open(input_files['gbk'], 'r') as gbk_file:
                        for line in gbk_file:
                            if 'transposase' in line:
                                summary.write(line)
                else:
                    summary.write("No Transposases\n")
                

                if len(open(input_files['amr']).readlines()) == 1:
                    summary.write("No AMR genes\n")
                else:
                    summary.write("AMR genes found\n")

                if len(open(input_files['vfdb']).readlines()) == 1:
                    summary.write("No virulence factor genes\n")
                else:
                    summary.write("Virulence factor genes found\n")

                if len(open(input_files['spacers']).readlines()) > 2:
                    summary.write("CRISPR spacers found\n")
                else:
                    summary.write("No CRISPR spacers found\n")
            else:
                summary.write("Genome includes multiple contigs, fragmented\n")
        else:
            summary.write("No contigs from the assembly were assigned viral, likely contigs too short in size\n")


def analyze_assembly(input_files, output_summary, params):
    file_content = Path(input_files['assembly']).read_text().splitlines()[-1]
    if 'Final assembly' in file_content or 'ALL DONE' in file_content:
        if Path(input_files['table']).exists():
            line_count = sum(1 for line in open(input_files['table']))
            if line_count > 1:
                copy_files(input_files, params)
                generate_summary(input_files, output_summary, params)
            else:
                with open(output_summary, 'w') as summary:
                    summary.write(f"Sample: {params['sample']}\n")
                    summary.write("Genome includes multiple contigs, fragmented\n")
        else:
            with open(output_summary, 'w') as summary:
                summary.write(f"Sample: {params['sample']}\n")
                summary.write("No contigs from the assembly were assigned viral, likely contigs too short in size\n")
    else:
        with open(output_summary, 'w') as summary:
            summary.write(f"Sample: {params['sample']}\n")
            summary.write("Failed during assembly\n")

# Replace input_files and output_params with the actual paths to your input/output files and parameters
input_files = {
        'assembly': snakemake.input.assembly,
        'table': snakemake.input.table,
        'genome': snakemake.input.genome,
        'gbk': snakemake.input.gbk,
        'plot': snakemake.input.plot,
        'amr': snakemake.input.amr,
        'vfdb': snakemake.input.vfdb,
        'spacers': snakemake.input.spacers,
        'taxa': snakemake.input.ph_taxa,
        'cdden': snakemake.input.cdden,
        'cds': snakemake.input.cds
}

output_summary = snakemake.output.summary

params = {
    'sample' : snakemake.params.sample,
    'genomes': snakemake.params.genomes,
    'gbks' : snakemake.params.gbks,
    'plots': snakemake.params.plots
}

# input_files = {
#     'assembly': '/home/nala0006/scratch/spae-paper/nanopore/Ecoli.out/PROCESSING/assembly/flye/Ecoli_17-sr/flye.log',
#     'table': '/home/nala0006/scratch/spae-paper/nanopore/Ecoli.out/PROCESSING/genome/Ecoli_17-sr/Ecoli_17-genome-candidates.csv',
#     'genome': '/home/nala0006/scratch/spae-paper/nanopore/Ecoli.out/PROCESSING/genome/Ecoli_17-sr/Ecoli_17_genome.fasta',
#     'gbk': '/home/nala0006/scratch/spae-paper/nanopore/Ecoli.out/PROCESSING/pharokka/Ecoli_17-sr/phynteny/phynteny.gbk',
#     'plot': '/home/nala0006/scratch/spae-paper/nanopore/Ecoli.out/PROCESSING/pharokka/Ecoli_17-sr/Ecoli_17_pharokka_plot.png',
#     'amr': '/home/nala0006/scratch/spae-paper/nanopore/Ecoli.out/PROCESSING/pharokka/Ecoli_17-sr/top_hits_card.tsv',
#     'vfdb': '/home/nala0006/scratch/spae-paper/nanopore/Ecoli.out/PROCESSING/pharokka/Ecoli_17-sr/top_hits_vfdb.tsv',
#     'spacers': '/home/nala0006/scratch/spae-paper/nanopore/Ecoli.out/PROCESSING/pharokka/Ecoli_17-sr/Ecoli_17_minced_spacers.txt'
# }

# output_summary = '/home/nala0006/example/Ecoli_17_summary.txt'

# params = {
#     'sample': 'Ecoli_17',
#     'dirout': Path('/home/nala0006/example/')
# }

analyze_assembly(input_files, output_summary, params)
