#!/usr/bin/env python3

import sys
import re
import gzip

def convert_chr_names(input_file, output_file):
    """
    Convert chromosome names from 'chr1' format to '1' format in a VCF file.
    Handles both gzipped (.gz) and uncompressed files.
    
    Args:
        input_file (str): Path to input VCF file (can be .vcf or .vcf.gz)
        output_file (str): Path to output VCF file (can be .vcf or .vcf.gz)
    """
    # Determine if files are gzipped based on file extension
    is_input_gzipped = input_file.endswith('.gz')
    is_output_gzipped = output_file.endswith('.gz')
    
    # Open input file with appropriate method
    infile = gzip.open(input_file, 'rt') if is_input_gzipped else open(input_file, 'r')
    
    # Open output file with appropriate method
    outfile = gzip.open(output_file, 'wt') if is_output_gzipped else open(output_file, 'w')
    
    try:
        for line in infile:
            # Preserve header lines
            if line.startswith('#'):
                outfile.write(line)
                continue
            
            # Convert chromosome names in data lines
            parts = line.split('\t')
            if parts[0].startswith('chr'):
                parts[0] = re.sub(r'^chr([0-9XYM])', r'\1', parts[0])
                line = '\t'.join(parts)
            
            outfile.write(line)
    finally:
        infile.close()
        outfile.close()

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: python vcf_chr_converter.py input.vcf[.gz] output.vcf[.gz]")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    try:
        convert_chr_names(input_file, output_file)
        print(f"Successfully converted chromosome names in {input_file}")
        print(f"Output written to {output_file}")
    except Exception as e:
        print(f"Error: {str(e)}")
        sys.exit(1)
