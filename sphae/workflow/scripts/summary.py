from pathlib import Path
import shutil
import pandas as pd
from Bio import SeqIO

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
                    summary.write(f"GC percent: {gc_percent}\n")
                    summary.write(f"Coding density: {coding_density}\n")

                found_gene = False
                gbk_records = SeqIO.parse(input_files['gbk'], 'genbank')
                for record in gbk_records:
                    for feature in record.features:
                        if feature.type == "CDS" and 'product' in feature.qualifiers and 'integra' in feature.qualifiers['product'][0].lower():
                            #print (feature)
                            found_gene = True
                            gene_id = feature.qualifiers['locus_tag'][0]
                            function = feature.qualifiers['function'][0]
                            product = feature.qualifiers['product'][0]
                            summary.write(f"\t{gene_id}: function=\"{function}\", product=\"{product}\"\n")

                if not found_gene:
                    if 'integra' not in open(input_files['gbk']).read().lower():
                        summary.write("No Integrases\n")
                    else:
                        summary.write("No Integrases\n")
                        summary.write("\t...but Phynteny predicted a few unknown function genes to have some similarity with integrase genes but with low confidence. Maybe a false positive or a novel integrase gene\n")

                    if 'recombinase' not in open(input_files['gbk']).read().lower():
                        summary.write("No recombinase\n")
                    else: 
                        summary.write("Recombinases found in genome\n")
                    if 'transposase' not in open(input_files['gbk']).read().lower():
                        summary.write("No transposase\n")
                    else: 
                        summary.write("Transposases found in genome\n")
                
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

"""
input_files = {
    'assembly': '/home/nala0006/scratch/wine_achromobacter_sarah_output/sphae.out/PROCESSING/assembly/flye/SQK-RBK114-24_barcode19-sr/flye.log',
    'table': '/home/nala0006/scratch/wine_achromobacter_sarah_output/sphae.out/PROCESSING/genome/SQK-RBK114-24_barcode19-sr/SQK-RBK114-24_barcode19-genome-candidates.csv',
    'genome': '/home/nala0006/scratch/wine_achromobacter_sarah_output/sphae.out/PROCESSING/genome/SQK-RBK114-24_barcode19-sr/SQK-RBK114-24_barcode19_genome.fasta',
    'gbk': '/home/nala0006/scratch/wine_achromobacter_sarah_output/sphae.out/PROCESSING/pharokka/SQK-RBK114-24_barcode19-sr/phynteny/phynteny.gbk',
    'plot': '/home/nala0006/scratch/wine_achromobacter_sarah_output/sphae.out/PROCESSING/pharokka/SQK-RBK114-24_barcode19-sr/SQK-RBK114-24_barcode19_pharokka_plot.png',
    'amr': '/home/nala0006/scratch/wine_achromobacter_sarah_output/sphae.out/PROCESSING/pharokka/SQK-RBK114-24_barcode19-sr/top_hits_card.tsv',
    'vfdb': '/home/nala0006/scratch/wine_achromobacter_sarah_output/sphae.out/PROCESSING/pharokka/SQK-RBK114-24_barcode19-sr/top_hits_vfdb.tsv',
    'taxa': '/home/nala0006/scratch/wine_achromobacter_sarah_output/sphae.out/PROCESSING/pharokka/SQK-RBK114-24_barcode19-sr/SQK-RBK114-24_barcode19_top_hits_mash_inphared.tsv',
    'spacers': '/home/nala0006/scratch/wine_achromobacter_sarah_output/sphae.out/PROCESSING/pharokka/SQK-RBK114-24_barcode19-sr/SQK-RBK114-24_barcode19_minced_spacers.txt',
    'cdden': '/home/nala0006/scratch/wine_achromobacter_sarah_output/sphae.out/PROCESSING/pharokka/SQK-RBK114-24_barcode19-sr/SQK-RBK114-24_barcode19_length_gc_cds_density.tsv',
    'cds': '/home/nala0006/scratch/wine_achromobacter_sarah_output/sphae.out/PROCESSING/pharokka/SQK-RBK114-24_barcode19-sr/SQK-RBK114-24_barcode19_cds_functions.tsv'
}

output_summary = '/home/nala0006/scratch/wine_achromobacter_sarah_output/sphae.out/RESULTS/SQK-RBK114-24_barcode19_summary.txt'

params = {
    'sample': 'SQK-RBK114-24_barcode19',
    'genomes': '/home/nala0006/scratch/wine_achromobacter_sarah_output/sphae.out/RESULTS/SQK-RBK114-24_barcode19-sr/SQK-RBK114-24_barcode19_genome.fasta',
    'gbks': '/home/nala0006/scratch/wine_achromobacter_sarah_output/sphae.out/RESULTS/SQK-RBK114-24_barcode19-sr/SQK-RBK114-24_barcode19.gbk',
    'plots': '/home/nala0006/scratch/wine_achromobacter_sarah_output/sphae.out/RESULTS/SQK-RBK114-24_barcode19_pharokka_plot.png'
}
"""
analyze_assembly(input_files, output_summary, params)
