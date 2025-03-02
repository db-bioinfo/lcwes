#!/usr/bin/env python3

def convert_chr_names(input_file, output_file):
    """
    Convert chromosome names from 1,2,3... to chr1,chr2,chr3...
    Handles regular chromosomes (1-22), X, Y, and M/MT
    """
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        # Process header line first
        header = infile.readline()
        outfile.write(header)
        
        # Process remaining lines
        for line in infile:
            parts = line.split('\t')
            if parts[0].isdigit():  # Regular chromosomes 1-22
                parts[0] = f'chr{parts[0]}'
            elif parts[0] in ['X', 'Y']:  # Sex chromosomes
                parts[0] = f'chr{parts[0]}'
            elif parts[0] in ['M', 'MT']:  # Mitochondrial
                parts[0] = 'chrM'
            
            # Write modified line
            outfile.write('\t'.join(parts))

if __name__ == '__main__':
    import sys
    
    if len(sys.argv) != 3:
        print("Usage: python script.py input_file output_file")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    convert_chr_names(input_file, output_file)
