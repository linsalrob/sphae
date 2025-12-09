import pandas as pd

def create_summary(pharokka_func, phold_func, pkl_func, tmp, summary_gbk):
    import pandas as pd

def create_summary(pharokka_func, phold_func, pkl_func, tmp, summary_gbk):
    pharokka_df = pd.read_csv(pharokka_func, sep='\t', header=None)
    phold_df    = pd.read_csv(phold_func,    sep='\t', header=None)
    pkl_df      = pd.read_csv(pkl_func,      sep='\t', header=None)

    # Find the maximum number of rows
    max_len = max(len(pharokka_df), len(phold_df), len(pkl_df))

    # Reindex all to the same length (pads with NaN)
    pharokka_df = pharokka_df.reindex(range(max_len))
    phold_df    = phold_df.reindex(range(max_len))
    pkl_df      = pkl_df.reindex(range(max_len))

    # Optional: save intermediate combined table
    tmp_df = pd.concat([pharokka_df, phold_df[2], pkl_df[2]], axis=1)
    tmp_df.to_csv(tmp, sep='\t', header=False, index=False)

    # Build final summary (NaN will appear automatically where missing)
    summary_df = pd.DataFrame({
        "contig name": pharokka_df[0],
        "protein ID":  pharokka_df[1],
        "pharokka":    pharokka_df[2],
        "phold":       phold_df[2],
        "phynteny":    pkl_df[2],
    })

    summary_df.columns = ["contig name", "protein ID", "pharokka", "phold", "phynteny"]
    # If you prefer literal "NA" instead of blank fields:
    summary_df = summary_df.fillna("NA")
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