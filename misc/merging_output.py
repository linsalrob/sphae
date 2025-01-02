import os
import glob
import pandas as pd
import argparse

def parse_summary(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    
    data = {
        "Sample": None,
        "Total length of reads after QC and subsampling": None,
        "Length": None,
        "Circular": None,
        "Graph connections": None,
        "Completeness": None,
        "Contamination": None,
        "Taxa Description": None,
        "Lowest Taxa classification": None,
        "Isolation host of described taxa": None,
        "Number of CDS": None,
        "Total number of CDS annotated as 'hypothetical protein'": None,
        "GC percent": None,
        "Coding density": None,
        "Gene search results": None,
        "Failed during assembly": None,
        "No contigs assigned viral": None
    }
    
    for line in lines:
        print(f"Processing line: {line.strip()}")
        if line.startswith("Sample:"):
            data["Sample"] = line.split(":")[1].strip()
        elif line.startswith("Total length of reads after QC and subsampling:"):
            data["Total length of reads after QC and subsampling"] = line.split(":")[1].strip()
        elif line.startswith("Length:"):
            data["Length"] = line.split(":")[1].strip()
        elif line.startswith("Circular:"):
            data["Circular"] = line.split(":")[1].strip()
        elif line.startswith("Graph connections:"):
            data["Graph connections"] = line.split(":")[1].strip()
        elif line.startswith("Completeness:"):
            data["Completeness"] = line.split(":")[1].strip()
        elif line.startswith("Contamination:"):
            data["Contamination"] = line.split(":")[1].strip()
        elif line.startswith("Taxa Description"):
            data["Taxa Description"] = line.split(":")[1].strip()
        elif line.startswith("Lowest Taxa classification:"):
            data["Lowest Taxa classification"] = line.split(":")[1].strip()
        elif line.startswith("Isolation host of described taxa:"):
            data["Isolation host of described taxa"] = line.split(":")[1].strip()
        elif line.startswith("Number of CDS:"):
            data["Number of CDS"] = line.split(":")[1].strip()
        elif line.startswith("Total number of CDS annotated as 'hypothetical protein':"):
            data["Total number of CDS annotated as 'hypothetical protein'"] = line.split(":")[1].strip()
        elif line.startswith("GC percent:"):
            data["GC percent"] = line.split(":")[1].strip()
        elif line.startswith("Coding density:"):
            data["Coding density"] = line.split(":")[1].strip()
        elif line.startswith("Gene search results:"):
            data["Gene search results"] = line.split(":")[1].strip()
        elif line.startswith("Failed during assembly"):
            data["Failed during assembly"] = "Yes"
        elif line.startswith("No contigs from the assembly were assigned viral, likely contigs too short in size"):
            data["No contigs assigned viral"] = "Yes"

    print(f"Parsed data: {data}")
    return data

def merge_summaries(directory, output_file):
    print(f"Directory: {directory}")
    summary_files = glob.glob(os.path.join(directory, "*_summary.txt"))
    print(f"Found summary files: {summary_files}")
    all_data = []

    for file_path in summary_files:
        print(f"Processing file: {file_path}")
        summary_data = parse_summary(file_path)
        all_data.append(summary_data)
    
    df = pd.DataFrame(all_data)
    print(f"DataFrame:\n{df}")
    df.to_csv(output_file, index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merge summary files into a single CSV file.")
    parser.add_argument("directory", help="Directory path containing summary files")
    parser.add_argument("output_file", help="Output CSV file name")
    args = parser.parse_args()

    merge_summaries(args.directory, args.output_file)