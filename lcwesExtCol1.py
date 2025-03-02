#!/usr/bin/env python3

import sys
import csv

def filter_fields(input_file, output_file):
    # Fields to keep (in order of appearance in original file)
    required_fields = [
        '#Chr', 'Start', 'End', 'Ref', 'Alt', 'avsnp151', 'Ref.Gene', 
        'Func.refGene', 'ExonicFunc.refGene',
        'clinvar: Clinvar ', # Note the extra space at the end
        ' InterVar: InterVar and Evidence ', # Note the spaces at both ends
        'Freq_gnomAD_genome_ALL', 'CADD_phred', 'SIFT_score',
        'dbscSNV_ADA_SCORE', 'dbscSNV_RF_SCORE', 'MetaSVM_score',
        'OMIM', 'Phenotype_MIM', 'OrphaNumber', 'Orpha', 'Otherinfo'
    ]

    try:
        with open(input_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
            # Read the header line to get column indices
            header = infile.readline().strip().split('\t')
            
            # Get indices of required fields
            field_indices = []
            for field in required_fields:
                try:
                    index = header.index(field)
                    field_indices.append(index)
                except ValueError:
                    print(f"Warning: Field '{field}' not found in input file", file=sys.stderr)
                    continue

            # Write header
            writer = csv.writer(outfile, delimiter='\t')
            writer.writerow([header[i] for i in field_indices])

            # Process data lines
            for line in infile:
                data = line.strip().split('\t')
                filtered_data = [data[i] for i in field_indices]
                writer.writerow(filtered_data)

        print(f"Processing complete. Output written to {output_file}", file=sys.stderr)

    except FileNotFoundError:
        print(f"Error: Input file {input_file} not found", file=sys.stderr)
        sys.exit(1)
    except PermissionError:
        print(f"Error: Permission denied when accessing files", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error: An unexpected error occurred: {str(e)}", file=sys.stderr)
        sys.exit(1)

def main():
    if len(sys.argv) != 3:
        print("Usage: python script.py input_file output_file", file=sys.stderr)
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    filter_fields(input_file, output_file)

if __name__ == "__main__":
    main()
