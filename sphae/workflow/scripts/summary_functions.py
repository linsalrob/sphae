import pandas as pd
import glob
import os

def extract_sample_name(file_path, suffix):
    return os.path.basename(file_path).replace(suffix, '').replace('_pharokka', '')

def create_summary(sample, pharokka_dir, phold_dir, pkl_dir, output_dir):
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

        # Check if any of the input files are empty
        if not (os.path.isfile(pharokka_func) and os.path.getsize(pharokka_func) > 0 
                and os.path.isfile(phold_func) and os.path.getsize(phold_func) > 0 
                and os.path.isfile(pkl_func) and os.path.getsize(pkl_func) > 0):
            # If any of the input files are empty, create an empty output file
            print(f"One or more input files for {sample_name} are empty. Creating empty output file.")
            output_file = os.path.join(output_dir, f"{sample_name}_summary.functions")
            open(output_file, 'w').close()  # Create an empty file
            continue

        # Load the data
        pharokka_df = pd.read_csv(pharokka_func, sep='\t', header=None)
        phold_df = pd.read_csv(phold_func, sep='\t', header=None)
        pkl_df = pd.read_csv(pkl_func, sep='\t', header=None)

        # Merge pharokka and phold on the second column (index 1)
        tmp_df = pd.merge(pharokka_df, phold_df, left_on=1, right_on=1, suffixes=('_pharokka', '_phold'))

        # Merge the intermediate results with pkl on the second column (index 1)
        final_df = pd.merge(tmp_df, pkl_df, left_on=1, right_on=1)

        # Select relevant columns
        summary_df = final_df.iloc[:, [0, 1, 2, 4, 6]].copy()
        summary_df.columns = ["contig name", "protein ID", "pharokka", "phold", "phynteny"]

        # Add sample name column using .assign
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
    sample=snakemake.params.ids,
)
