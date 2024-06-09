from pathlib import Path
import shutil
import pandas as pd
from Bio import SeqIO
import os
import glob

def copy_files(input_files, params):
    shutil.copy(input_files['genome'], params['genomes'])
    shutil.copy(input_files['gbk'], params['gbks'])
    shutil.copy(input_files['plot'], params['plots'])

def copy_multiple_files(params):
    # Ensure the destination directory exists
    os.makedirs(outdir, exist_ok=True)
    
    # Use glob to find all .fasta files in the specified pattern
    source_files = glob.glob(f"{params['annot']}/{params['sample']}-*-genomes/*.fasta")
    # Copy each .fasta file to the destination directory
    for file in source_files:
        print(f"Copying {file} to {outdir}")
        shutil.copy(file, outdir)
    
    # Find and copy .gbk files with renaming
    gbk_files = glob.glob(f"{params['annot']}/phynteny-*/{params['sample']}_*_phynteny/*.gbk")
    for i, gbk in enumerate(gbk_files, start=1):
        samplenames = f"{params['sample']}_{i}"
        new_gbk_path = os.path.join(outdir, f"{samplenames}.gbk")
        print(f"Copying {gbk} to {new_gbk_path}")
        shutil.copy(gbk, new_gbk_path)
        i=i+1

    # Find and copy _plot.png files with renaming
    plot_files = glob.glob(f"{params['annot']}/phynteny-*/{params['sample']}_*_phynteny/*_plot.png")
    for i, plt in enumerate(plot_files, start=1):
        samplenames = f"{params['sample']}_{i}"
        new_png_path = os.path.join(outdir, f"{samplenames}.png")
        print(f"Copying {plt} to {new_png_path}")
        shutil.copy(plt, new_png_path)
        i=i+1
    
def write_single_genome_summary(input_files, summary):
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
            
            #AMR genes 
            if (len(open(input_files['amr']).readlines()) == 1) and (len(open(input_files['card']).readlines()) == 0):
                summary.write("No AMR genes found\n")
            else:
                summary.write("AMR genes found\n")                

            #Virulence genes 
            if (len(open(input_files['vfdb']).readlines()) == 1) and (len(open(input_files['vfdb_phold']).readlines()) == 0):
                summary.write("No virulence factor genes\n")
            else:
                summary.write("Virulence genes found\n")

            #CRISPR spacers 
            if (len(open(input_files['spacers']).readlines()) == 0) and (len(open(input_files['acr']).readlines()) == 0):
                summary.write("No CRISPR spacers found\n")
            else:
                summary.write("CRISPR spacers found\n")
                
            #Defense finder genes  
            if (len(open(input_files['defense']).readlines()) == 0):
                summary.write("No Defense genes found\n")
            else:
                summary.write("Defense genes found\n")   

def write_multiple_genome_summary(input_files, summary):
    summary.write("Multiple phages assembled from this sample\n")
    summary.write("Their characteristics are:\n")
    summary.write("\n\n")

    with open(input_files['table'], 'r') as table_file:
        lines = table_file.readlines()
        if len(lines) > 2:
            for i, line in enumerate(lines[1:], start=1):  # Skip header and enumerate starting from 1
                fields = line.split(',')
                samplenames = f"{params['sample']}_{i}"
                summary.write(f"Sample name: {samplenames}\n")
                summary.write(f"Length: {fields[2]}\n")
                if '+' in fields[8]:
                    summary.write("Circular: True\t")
                else:
                    summary.write("Circular: False\n")

                if 'DTR' in fields[21]:
                    summary.write("DTR found\n")

                summary.write(f"Graph connections: {fields[4]}\n")
                summary.write(f"Completeness: {fields[20]}\n")
                summary.write(f"Contamination: {fields[22]}\n")
                
                taxa_pattern = f"{annot}/pharokka-pr/{samplenames}_pharokka/{samplenames}_top_hits_mash_inphared.tsv"
                for taxa_file in glob.glob(taxa_pattern):
                    tax = pd.read_csv(taxa_file, sep='\t')
                    summary.write("Taxa name (Matching hashes):\t")
                    for index, row in tax.iterrows():
                        summary.write(f"{row['Description']}\t{row['mash_matching_hashes']}\n")

                cds_pattern = f"{annot}/pharokka-pr/{samplenames}_pharokka/{samplenames}_cds_functions.tsv"
                for cds_file in glob.glob(cds_pattern):
                    cds_df=pd.read_csv(cds_file, sep='\t')
                    cds_data = cds_df[cds_df['Description'] == 'CDS']
                    count_value = cds_data['Count'].values[0]
                    summary.write(f"Number of CDS: {count_value}\n")
                                 
                cds_dens_pattern = f"{annot}/pharokka-pr/{samplenames}_pharokka/{samplenames}_length_gc_cds_density.tsv"
                for cdsdens_file in glob.glob(cds_dens_pattern):
                    cdn=pd.read_csv(cdsdens_file, sep='\t')
                    gc_percent = cdn['gc_perc'].values[0]
                    coding_density = cdn['cds_coding_density'].values[0]
                    summary.write(f"GC percent: {gc_percent}\n")
                    summary.write(f"Coding density: {coding_density}\n")

                found_gene = False
                gbk_pattern = f"{annot}/phynteny-pr/{samplenames}_phynteny/phynteny.gbk"
                for gbk_file in glob.glob(gbk_pattern):
                    gbk_records = SeqIO.parse(gbk_file, 'genbank')
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
                        if 'integra' not in open(gbk_file).read().lower():
                            summary.write("No Integrases\n")
                        else:
                            summary.write("No Integrases\n")
                            summary.write("\t...but Phynteny predicted a few unknown function genes to have some similarity with integrase genes but with low confidence. Maybe a false positive or a novel integrase gene\n")

                        if 'recombinase' not in open(gbk_file).read().lower():
                            summary.write("No recombinase\n")
                        else: 
                            summary.write("Recombinases found in genome\n")
                        
                        if 'transposase' not in open(gbk_file).read().lower():
                            summary.write("No transposase\n")
                        else: 
                            summary.write("Transposases found in genome\n")
                        
                #AMR genes 
                amr_pattern = f"{annot}/pharokka-pr/{samplenames}_pharokka/top_hits_card.tsv"
                card_pattern = f"{annot}/phold-pr/{samplenames}_phold/sub_db_tophits/card_cds_predictions.tsv"
                amr_found = any(glob.glob(amr_pattern)) or any(glob.glob(card_pattern))
                summary.write("AMR genes found\n" if amr_found else "No AMR genes found\n")
               

                # Virulence genes
                vfdb_pattern = f"{annot}/pharokka-pr/*/{samplenames}_pharokka/vfdb_annotations.tsv"
                vfdb_phold_pattern = f"{annot}/pharokka-pr/*/{samplenames}_pharokka/phage_vfdb_annotations.tsv"
                virulence_found = any(glob.glob(vfdb_pattern)) or any(glob.glob(vfdb_phold_pattern))
                summary.write("Virulence genes found\n" if virulence_found else "No virulence factor genes\n")

                # CRISPR spacers
                spacers_pattern = f"{annot}/pharokka-pr/*/{samplenames}_pharokka/crispr_summary.tsv"
                acr_pattern = f"{annot}/pharokka-pr/*/{samplenames}_pharokka/acr_annotations.tsv"
                spacers_found = any(glob.glob(spacers_pattern)) or any(glob.glob(acr_pattern))
                summary.write("CRISPR spacers found\n" if spacers_found else "No CRISPR spacers found\n")

                # Defense genes
                defense_pattern = f"{annot}/pharokka-pr/*/{samplenames}_pharokka/defense_summary.tsv"
                defense_found = any(glob.glob(defense_pattern))
                summary.write("Defense genes found\n" if defense_found else "No Defense genes found\n")
                summary.write("\n\n")

def generate_summary(input_files, output_summary, params):
    with open(output_summary, 'w') as summary:
        summary.write(f"Sample: {params['sample']}\n")
        if Path(input_files['table']).exists():
            line_count = sum(1 for line in open(input_files['table']))
            if line_count == 2:
                write_single_genome_summary(input_files, summary)
            elif line_count > 2:
                write_multiple_genome_summary(input_files, summary)
        else:
            summary.write("No contigs from the assembly were assigned viral, likely contigs too short in size\n")

def analyze_assembly(input_files, output_summary, params):
    file_content = Path(input_files['assembly']).read_text().splitlines()[-1]
    if 'Final assembly' in file_content or 'ALL DONE' in file_content:
        if Path(input_files['table']).exists():
            line_count = sum(1 for line in open(input_files['table']))
            if line_count == 2:
                copy_files(input_files, params)
                generate_summary(input_files, output_summary, params)
            elif line_count > 2:
                copy_multiple_files(params)
                generate_summary(input_files, output_summary, params)
            elif line_count == 1 :
                with open(output_summary, 'w') as summary:
                    summary.write(f"Sample: {params['sample']}\n")
                    summary.write("No contigs from the assembly were assigned viral, likely contigs too short in size\n")
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
        'taxa': snakemake.input.ph_taxa,
        'cdden': snakemake.input.cdden,
        'cds': snakemake.input.cds,
        'amr': snakemake.input.amr,
        'vfdb': snakemake.input.vfdb,
        'spacers': snakemake.input.spacers,
        'acr': snakemake.input.acr,
        'card': snakemake.input.card,
        'defense': snakemake.input.defense,
        'vfdb_phold': snakemake.input.vfdb_phold,
}

output_summary = snakemake.output.summary
annot= snakemake.params.annot
outdir=snakemake.params.outdir
params = {
    'sample' : snakemake.params.ID,
    'genomes': snakemake.params.genomes,
    'gbks' : snakemake.params.gbks,
    'plots': snakemake.params.plots,
    'annot': snakemake.params.annot
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
