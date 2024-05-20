import pandas as pd

def create_summary(pharokka_func, phold_func, pkl_func, tmp, summary_gbk):
    # Load the data
    pharokka_df = pd.read_csv(pharokka_func, sep='\t', header=None)
    phold_df = pd.read_csv(phold_func, sep='\t', header=None)
    pkl_df = pd.read_csv(pkl_func, sep='\t', header=None)

    # Merge pharokka and phold on the second column (index 1)
    tmp_df = pd.merge(pharokka_df, phold_df, left_on=1, right_on=1, suffixes=('_pharokka', '_phold'))

    # Save intermediate results
    tmp_df.to_csv(tmp, sep='\t', header=False, index=False)

    # Merge the intermediate results with pkl on the second column (index 1)
    final_df = pd.merge(tmp_df, pkl_df, left_on=1, right_on=1)

    # Select relevant columns
    summary_df = final_df.iloc[:, [0, 1, 2, 4, 6]]
    summary_df.columns = ["contig name", "protein ID", "pharokka", "phold", "phynteny"]

    # Save the final summary
    summary_df.to_csv(summary_gbk, sep='\t', header=True, index=False)

# Example usage
create_summary(
    pharokka_func=snakemake.input.pharokka_func,
    phold_func=snakemake.input.phold_func,
    pkl_func=snakemake.input.pkl_func,
    tmp=snakemake.params.tmp,
    summary_gbk=snakemake.output.summary_gbk
)
