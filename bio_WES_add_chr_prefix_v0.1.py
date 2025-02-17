#!/usr/bin/env python3

import sys
import csv

def convert_chromosome_names(input_file, output_file):
    try:
        with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
            # Read the header line
            header = infile.readline()
            outfile.write(header)
            
            # Process the rest of the file
            for line in infile:
                # Split the line by tabs
                fields = line.split('\t')
                
                # Check if the first field is a chromosome number
                if fields[0].isdigit():
                    # Add 'chr' prefix
                    fields[0] = 'chr' + fields[0]
                elif fields[0].upper() in ['X', 'Y', 'M', 'MT']:
                    # Handle X, Y, and mitochondrial chromosomes
                    fields[0] = 'chr' + fields[0].upper()
                
                # Join the fields back together and write to output
                outfile.write('\t'.join(fields))
                
        print(f"Successfully converted {input_file} to {output_file}")
        
    except FileNotFoundError:
        print(f"Error: Could not find the input file {input_file}")
        sys.exit(1)
    except PermissionError:
        print(f"Error: Permission denied when trying to write to {output_file}")
        sys.exit(1)
    except Exception as e:
        print(f"An unexpected error occurred: {str(e)}")
        sys.exit(1)

def main():
    if len(sys.argv) != 3:
        print("Usage: python script.py input_file output_file")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    convert_chromosome_names(input_file, output_file)

if __name__ == "__main__":
    main()
