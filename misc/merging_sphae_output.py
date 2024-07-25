#!/usr/bin/env python3

import os
import pandas as pd
import argparse

# Define function to parse a summary file
def parse_summary_file(file_path):
    entries = []
    entry = {}
    multiple_phages = False
    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if 'Sample:' in line:
                # Extract the sample name
                key, value = line.split(':', 1)
                entry['sample'] = value.strip()
            elif "No contigs from the assembly were assigned viral" in line:
                # Custom handling for the special case
                entry['length'] = "No contigs from the assembly were assigned viral, likely contigs too short in size"
            elif 'Multiple phages assembled from this sample' in line:
                # Start of multiple phage entries
                multiple_phages = True
                if entry:
                    entries.append(entry)
            elif line.startswith('Sample name:') and multiple_phages:
                # New phage entry within multiple phages
                if entry:
                    entries.append(entry)
                entry = {'sample': entry['sample']}
                entry['sample_name'] = line.split(':', 1)[1].strip()  # Ensure line contains ':' before splitting
            elif ':' in line:
                # Standard key-value pair
                key, value = line.split(':', 1)
                key = key.strip().lower().replace(' ', '_').replace('(', '').replace(')', '')
                value = value.strip()
                entry[key] = value
            else:
                # Special handling for lines without a colon
                line_cleaned = line.strip().lower().replace(' ', '_').replace('(', '').replace(')', '')
                entry[line_cleaned] = True

        # Append the last entry if any
        if entry:
            entries.append(entry)

    return entries

def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='Parse summary files and convert to a DataFrame.')
    parser.add_argument('base_dir', type=str, help='Base directory containing summary files.')
    parser.add_argument('output_file', type=str, help='Output file name for the resulting DataFrame.')

    args = parser.parse_args()

    # List to hold parsed data
    summary_data = []

    # Traverse directories and read summary files
    for dir_name in os.listdir(args.base_dir):
        dir_path = os.path.join(args.base_dir, dir_name)
        if os.path.isdir(dir_path):  # Check if it is a directory
            # Read files in the directory
            for file_name in os.listdir(dir_path):
                if file_name.endswith('_summary.txt'):
                    file_path = os.path.join(dir_path, file_name)
                    # Parse the summary file and extend summary_data with parsed entries
                    summary_data.extend(parse_summary_file(file_path))

    # Convert to DataFrame
    df = pd.DataFrame(summary_data)

    # Write DataFrame to a TSV file
    df.to_csv(args.output_file, sep='\t', index=False)

    print(f"DataFrame has been successfully written to '{args.output_file}'.")

if __name__ == '__main__':
    main()
