#!/usr/bin/env python3

import argparse
import gzip
import re
import sys

def parse_info_field(info_str):
    """Extract required fields from the INFO field"""
    info_dict = {}
    items = info_str.split(';')
    
    for item in items:
        if '=' in item:
            key, value = item.split('=', 1)
            info_dict[key] = value
        else:
            info_dict[item] = True
    
    return info_dict

def extract_vcf_data(vcf_file, output_file):
    """Extract required fields from VCF file and write to TSV file"""
    
    # Open output file
    with open(output_file, 'w') as out_fh:
        # Write header with VAF (even though we're extracting AF)
        out_fh.write("Chr\tStart\tEnd\tRef\tAlt\tGene\tFilter\tCLNHGVS\tVAF\tRS\n")
        
        # Open VCF file (gzipped)
        opener = gzip.open if vcf_file.endswith('.gz') else open
        with opener(vcf_file, 'rt') as in_fh:
            for line in in_fh:
                # Skip header lines
                if line.startswith('#'):
                    continue
                
                # Parse VCF line
                fields = line.strip().split('\t')
                if len(fields) < 8:
                    sys.stderr.write(f"Warning: Skipping malformed line: {line.strip()}\n")
                    continue
                
                # Extract basic fields
                chrom = fields[0]
                pos = fields[1]
                rs_id = fields[2]  # Use RS ID from VCF
                ref = fields[3]
                alt = fields[4]
                filter_val = fields[6]  # FILTER column
                
                # Calculate end position (start + length of reference - 1)
                try:
                    end = str(int(pos) + len(ref) - 1)
                except ValueError:
                    sys.stderr.write(f"Warning: Invalid position value: {pos}\n")
                    continue
                
                # Parse INFO field
                info_dict = parse_info_field(fields[7])
                
                # Extract CLNHGVS
                clnhgvs = info_dict.get('CLNHGVS', 'NA')
                
                # Extract AF from INFO field first (HaplotypeCaller puts it there)
                af = info_dict.get('AF', 'NA')
                
                # Extract RS if not already present in the ID column
                rs = 'NA'
                if rs_id != '.' and rs_id.startswith('rs'):
                    rs = rs_id
                elif 'RS' in info_dict:
                    rs = info_dict['RS']
                
                # Extract gene name
                gene = 'NA'
                # First check GENEINFO which is present in your VCF format
                if 'GENEINFO' in info_dict and ':' in info_dict['GENEINFO']:
                    gene = info_dict['GENEINFO'].split(':', 1)[0]
                # Then try ANN field which contains detailed annotations
                elif 'ANN' in info_dict:
                    ann_parts = info_dict['ANN'].split('|')
                    if len(ann_parts) > 3:  # Ensure we have enough parts
                        gene = ann_parts[3]  # Gene name is typically in 4th position
                    # If gene is still NA, try to extract it from the first annotation field
                    if gene == 'NA' and '|' in info_dict['ANN']:
                        first_ann = info_dict['ANN'].split(',')[0]
                        ann_fields = first_ann.split('|')
                        if len(ann_fields) > 3:
                            gene = ann_fields[3]
                # Finally try other common gene fields
                else:
                    for gene_field in ['GENE', 'Gene', 'gene', 'SYMBOL']:
                        if gene_field in info_dict:
                            gene = info_dict[gene_field]
                            break
                
                # If AF is still NA, check FORMAT and sample fields
                if af == 'NA' and len(fields) > 9:
                    format_fields = fields[8].split(':')
                    sample_fields = fields[9].split(':')
                    
                    # Look for AF field in FORMAT
                    if 'AF' in format_fields:
                        af_index = format_fields.index('AF')
                        if len(sample_fields) > af_index:
                            af = sample_fields[af_index]
                
                # Write output line
                out_line = f"{chrom}\t{pos}\t{end}\t{ref}\t{alt}\t{gene}\t{filter_val}\t{clnhgvs}\t{af}\t{rs}\n"
                out_fh.write(out_line)

def main():
    parser = argparse.ArgumentParser(description='Extract specific fields from VCF file')
    parser.add_argument('input_file', help='Input VCF file (can be gzipped)')
    parser.add_argument('output_file', help='Output TSV file')
    args = parser.parse_args()
    
    extract_vcf_data(args.input_file, args.output_file)
    print(f"Extraction complete. Results written to {args.output_file}")

if __name__ == "__main__":
    main()
