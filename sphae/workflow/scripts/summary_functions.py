import pandas as pd
import glob
import os

def extract_sample_name(file_path, suffix):
    return os.path.basename(file_path).replace(suffix, '').replace('_pharokka', '')

def create_summary(sample, pharokka_dir, phold_dir, pkl_dir, tmp, output_dir):
    # Find all *.functions files in the pharokka directory (assuming others follow the same structure)
    pharokka_funcs = glob.glob(os.path.join(pharokka_dir, f"{sample}_*_pharokka", '*_pharokka.functions'))

    # Extract sample names
    sample_names = [extract_sample_name(f, '_pharokka.functions') for f in pharokka_funcs]
    print("Sample names:", sample_names) 

    # Process each sample individually
    for sample_name in sample_names:
       # Define file paths
        pharokka_func = os.path.join(pharokka_dir, f"{sample_name}_pharokka", f"{sample_name}_pharokka.functions")
        phold_func = os.path.join(phold_dir, f"{sample_name}_phold", f"{sample_name}_phold.functions")
        pkl_func = os.path.join(pkl_dir, f"{sample_name}_phynteny", "phynteny.functions")
        
        # Load the data
        pharokka_df = pd.read_csv(pharokka_func, sep='\t', header=None, names=['contig', 'protein_id', 'pharokka'])
        phold_df = pd.read_csv(phold_func, sep='\t', header=None, names=['contig', 'protein_id', 'phold'])
        pkl_df = pd.read_csv(pkl_func, sep='\t', header=None, names=['contig', 'protein_id', 'phynteny'])

        # Find the maximum number of rows
        max_len = max(len(pharokka_df), len(phold_df), len(pkl_df))

        # Reindex all to the same length (pads with NaN)
        pharokka_df = pharokka_df.reindex(range(max_len))
        phold_df    = phold_df.reindex(range(max_len))
        pkl_df      = pkl_df.reindex(range(max_len))

        # Optional: save intermediate combined table
        tmp_df = pd.concat([pharokka_df, phold_df['phold'], pkl_df['phynteny']], axis=1)
        tmp_file = os.path.join(tmp, f"{sample_name}_tmp")
        tmp_df.to_csv(tmp, sep='\t', header=False, index=False)

        # Build final summary (NaN will appear automatically where missing)
        summary_df = pd.DataFrame({
            "contig name": pharokka_df['contig'],
            "protein ID":  pharokka_df['protein_id'],
            "pharokka":    pharokka_df['pharokka'],
            "phold":       phold_df['phold'],
            "phynteny":    pkl_df['phynteny'],
        })
        # If you prefer literal "NA" instead of blank fields:
        summary_df = summary_df.fillna("NA")
        # Save the final summary
        summary_df = summary_df.assign(sample_name=sample_name)
        # Save the final summary for each sample
        output_file = os.path.join(output_dir, f"{sample_name}_summary.functions")
        summary_df.to_csv(output_file, sep='\t', header=True, index=False)
        print(f"Saved summary for {sample_name} to {output_file}")

# Example usage
create_summary(
    pharokka_dir=snakemake.params.pharokka_func,
    phold_dir=snakemake.params.phold_func,
    pkl_dir=snakemake.params.pkl_func,
    output_dir=snakemake.params.output,
    tmp=snakemake.params.tmp,
    sample=snakemake.params.ids,
)
