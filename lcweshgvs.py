#!/usr/bin/env python3
import csv
import sys
import argparse

def main():
    # Set up command-line argument parsing
    parser = argparse.ArgumentParser(description='Enrich a variant file with HGVS and Type information from another file.')
    parser.add_argument('file1', help='Main variant file to be enriched')
    parser.add_argument('file2', help='File with HGVS and Type information')
    parser.add_argument('output', help='Output enriched file path')
    
    # Parse arguments
    args = parser.parse_args()
    
    # Define input and output files from arguments
    file1_path = args.file1  # Main file to be enriched
    file2_path = args.file2  # File with HGVS and Type information
    output_path = args.output  # Output file
    
    # Create a dictionary to store variants from file2
    # Key will be a tuple of (chr, pos, ref, alt)
    # Value will be a tuple of (HGVS, Type)
    file2_variants = {}
    
    # Read file2 and populate the dictionary
    print("Reading file2.tsv for HGVS and Type information...")
    with open(file2_path, 'r') as f2:
        reader = csv.reader(f2, delimiter='\t')
        for row in reader:
            # Skip any header row if it exists
            if row[0] == 'chr' or row[0].startswith('#'):
                continue
                
            # Extract the necessary fields
            chrom = row[0]
            pos = row[1]  # Start position
            ref = row[3]  # Reference allele
            alt = row[4]  # Alternative allele
            hgvs = row[5]  # HGVS (column 6)
            variant_type = row[6]  # Type (column 7)
            
            # Create the key for the dictionary
            key = (chrom, pos, ref, alt)
            
            # Store the HGVS and Type for this variant
            file2_variants[key] = (hgvs, variant_type)
    
    print(f"Loaded {len(file2_variants)} variants from file2.tsv")
    
    # Process file1 and create the enriched output
    print("Processing file1.tsv and adding HGVS and Type columns...")
    enriched_count = 0
    total_count = 0
    
    with open(file1_path, 'r') as f1, open(output_path, 'w', newline='') as f_out:
        reader = csv.reader(f1, delimiter='\t')
        writer = csv.writer(f_out, delimiter='\t')
        
        # Process each row in file1
        for row_num, row in enumerate(reader):
            total_count += 1
            
            # Handle the header row
            if row_num == 0:
                # Add new column headers
                row.extend(["HGVS", "Type"])
                writer.writerow(row)
                continue
            
            # Extract the necessary fields for matching
            chrom = row[0]
            pos = row[1]  # Start position
            ref = row[3]  # Reference allele
            alt = row[4]  # Alternative allele
            
            # Create the key for lookup
            key = (chrom, pos, ref, alt)
            
            # Look up the variant in our dictionary
            if key in file2_variants:
                # Get the HGVS and Type
                hgvs, variant_type = file2_variants[key]
                
                # Add the new columns
                row.extend([hgvs, variant_type])
                enriched_count += 1
            else:
                # If variant not found in file2, add empty values
                row.extend(["", ""])
            
            # Write the modified row to the output file
            writer.writerow(row)
    
    print(f"Processed {total_count} variants from file1.tsv")
    print(f"Enriched {enriched_count} variants with HGVS and Type information")
    print(f"Output written to {output_path}")

if __name__ == "__main__":
    # Check if arguments provided
    if len(sys.argv) == 1:
        print("Usage: python script.py file1.tsv file2.tsv output.tsv")
        print("Run with -h for more information")
        sys.exit(1)
    main()
