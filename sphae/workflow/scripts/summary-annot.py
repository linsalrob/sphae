from pathlib import Path
import shutil
import pandas as pd
from Bio import SeqIO
import glob, os

def copy_files(input_files, params):
    shutil.copy(input_files['genome'], params['genomes'])
    shutil.copy(input_files['gbk'], params['gbks'])

    os.makedirs(params['plots'], exist_ok=True)

    plot_files = glob.glob(os.path.join(input_files['plot'], "*.svg")) + \
                 glob.glob(os.path.join(input_files['plot'], "*.png"))

    if len(plot_files) == 0:
        print("Warning: No plot files found")
    else:
        for f in plot_files:
            shutil.copy(f, params['plots'])

def count_hypothetical_proteins(params):
    count = 0
    for record in SeqIO.parse(params['gbks'], "genbank"):
        for feature in record.features:
            if feature.type == "CDS":
                if "product" in feature.qualifiers:
                    if any("hypothetical protein" in p.lower() for p in feature.qualifiers["product"]):
                        count += 1
    return count

def is_protein(seq):
    return any(c not in "ACGTUNacgtun" for c in seq if c.isalpha())

def generate_summary(input_files, output_summary, params):
    copy_files(input_files, params)

    with open(output_summary, 'w') as summary:
        summary.write(f"Sample: {params['sample']}\n")

        # --- Sequence type detection (safe + fast) ---
        records = SeqIO.parse(input_files['genome'], "fasta")
        first_record = next(records, None)

        if first_record is None:
            summary.write("Empty genome file\n")
            return

        is_prot = is_protein(str(first_record.seq))
        summary.write(f"Sequence type: {'protein' if is_prot else 'nucleotide'}\n")

        # --- Contig count ---
        if not is_prot:
            contigs = sum(1 for _ in SeqIO.parse(input_files['genome'], "fasta"))
            summary.write(f"Number of contigs: {contigs}\n")

        # --- Taxa ---
        tax = pd.read_csv(input_files['taxa'], sep='\t')
        summary.write("Taxa name (Matching hashes):\n")
        for _, row in tax.iterrows():
            summary.write(f"{row['Description']}\t{row['mash_matching_hashes']}\n")

        # --- CDS count ---
        cds_df = pd.read_csv(input_files['cds'], sep='\t')
        cds_data = cds_df[cds_df['Description'] == 'CDS']
        if not cds_data.empty:
            summary.write(f"Number of CDS: {cds_data['Count'].values[0]}\n")

        # --- Hypothetical proteins ---
        hypothetical_protein_count = count_hypothetical_proteins(params)
        summary.write(f"Total number of CDS annotated as 'hypothetical protein': {hypothetical_protein_count}\n")

        # --- GC + coding density ---
        cdn = pd.read_csv(input_files['cdden'], sep='\t')
        summary.write(f"GC percent: {cdn['gc_perc'].values[0]}\n")
        summary.write(f"Coding density: {cdn['cds_coding_density'].values[0]}\n")

        # --- Read GBK once ---
        gbk_text = open(input_files['gbk']).read().lower()

        # --- Integrases ---
        found_gene = False
        for record in SeqIO.parse(input_files['gbk'], 'genbank'):
            for feature in record.features:
                if feature.type == "CDS" and 'product' in feature.qualifiers:
                    product = feature.qualifiers['product'][0].lower()
                    if 'integra' in product:
                        found_gene = True
                        gene_id = feature.qualifiers.get('locus_tag', ['unknown'])[0]
                        function = feature.qualifiers.get('function', ['unknown'])[0]
                        summary.write(f"\t{gene_id}: function=\"{function}\", product=\"{product}\"\n")

        if not found_gene:
            summary.write("No Integrases\n")
            if 'integra' in gbk_text:
                summary.write("\t...but possible low-confidence hits detected\n")

        # --- Other mobile elements ---
        summary.write("Recombinases found in genome\n" if 'recombinase' in gbk_text else "No recombinase\n")
        summary.write("Transposases found in genome\n" if 'transposase' in gbk_text else "No transposase\n")

        # --- AMR ---
        if (len(open(input_files['amr']).readlines()) == 1) and (len(open(input_files['card']).readlines()) == 0):
            summary.write("No AMR genes found\n")
        else:
            summary.write("AMR genes found\n")

        # --- Virulence ---
        if (len(open(input_files['vfdb']).readlines()) == 1) and (len(open(input_files['vfdb_phold']).readlines()) == 0):
            summary.write("No virulence factor genes\n")
        else:
            summary.write("Virulence genes found\n")

        # --- CRISPR ---
        if (len(open(input_files['spacers']).readlines()) == 0) and (len(open(input_files['acr']).readlines()) == 0):
            summary.write("No CRISPR spacers found\n")
        else:
            summary.write("CRISPR spacers found\n")

        # --- Defense ---
        if len(open(input_files['defense']).readlines()) == 0:
            summary.write("No Defense genes found\n")
        else:
            summary.write("Defense genes found\n")


# --- Snakemake integration ---
input_files = {
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

params = {
    'sample': snakemake.params.sample,
    'genomes': snakemake.params.genomes,
    'gbks': snakemake.params.gbks,
    'plots': snakemake.params.plots,
    #'input_type': snakemake.params.input_type
}

generate_summary(input_files, output_summary, params)