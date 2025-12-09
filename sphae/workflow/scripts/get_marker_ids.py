#!/usr/bin/env python3

import os
import glob

def parse_attrs(attr_str: str) -> dict:
    """Parse GFF3 9th column into a dict: key=value;key2=value2;..."""
    d = {}
    for item in attr_str.strip().split(";"):
        if not item:
            continue
        if "=" in item:
            k, v = item.split("=", 1)
            d[k] = v
    return d

def find_gff_files(pattern_or_dir: str):
    """Return a list of GFF3 files from either a glob pattern or a directory."""
    if os.path.isdir(pattern_or_dir):
        # walk directory and grab any *.gff3
        gffs = []
        for root, dirs, files in os.walk(pattern_or_dir):
            for f in files:
                if f.endswith(".gff3"):
                    gffs.append(os.path.join(root, f))
        return gffs
    else:
        # treat as glob pattern like 'phynteny-pr/*/*.gff3'
        return [p for p in glob.glob(pattern_or_dir) if os.path.isfile(p)]


# Snakemake params
pattern = snakemake.params.folder      # can be folder OR glob like 'phynteny-pr/*/*.gff3'
interest = snakemake.params.interest.lower()
out_fasta = snakemake.output[0] if isinstance(snakemake.output, list) else snakemake.output

gff_files = find_gff_files(pattern)

print("[DEBUG] GFF files I see:")
for p in gff_files:
    print("  -", p)

n_features = 0
n_hits = 0

with open(out_fasta, "w") as out:
    for path in gff_files:
        with open(path) as fh:
            for line in fh:
                if line.startswith("#") or "\t" not in line:
                    continue
                cols = line.rstrip("\n").split("\t")
                if len(cols) < 9:
                    continue

                attrs = parse_attrs(cols[8])
                product = attrs.get("product", "").lower()
                if interest not in product:
                    continue

                _id = attrs.get("ID") or attrs.get("protein_id")
                seq = attrs.get("translation")
                if not _id or not seq:
                    continue

                n_features += 1
                n_hits += 1
                out.write(f">{_id}\n{seq}\n")

print(
    f"Scanned {n_features} matching features in {len(gff_files)} GFF files; "
    f"wrote {n_hits} sequences to {out_fasta}"
)
