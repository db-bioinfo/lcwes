#!/usr/bin/env python3
"""
Script to add VCF annotations to a TSV file.
Matches variants between the TSV and VCF files (both GATK and FreeBayes)
and adds specified fields from the VCF as new columns in the TSV.
"""

import pandas as pd
import gzip
import re
import argparse
import sys
import os

def extract_info_field(info, field):
    """Extract a specific field from the INFO column."""
    match = re.search(f"{field}=([^;]+)", info)
    return match.group(1) if match else ""

def extract_dp(format_col, sample_col):
    """Extract DP value from FORMAT and sample columns."""
    format_fields = format_col.split(':')
    sample_values = sample_col.split(':')
    
    if 'DP' in format_fields:
        dp_index = format_fields.index('DP')
        if dp_index < len(sample_values):
            return sample_values[dp_index]
    return ""

def extract_snpeff_info(ann_field):
    """Extract effect and location from SnpEff ANN field."""
    if not ann_field:
        return "", ""
    
    annotations = ann_field.split(',')
    parts = annotations[0].split('|')  # Use the first annotation
    
    effect = parts[1] if len(parts) > 1 else ""
    
    # Extract location information (exon/intron number)
    location = ""
    if len(parts) > 8:
        rank_total = parts[8]
        if rank_total and '/' in rank_total:
            feature_type = "Exon"  # Default to Exon
            if "intron_variant" in effect:
                feature_type = "Intron"
            location = f"{feature_type} {rank_total}"
    
    return effect, location

def parse_vcf(vcf_file):
    """Parse VCF file and extract relevant information for each variant."""
    variants = {}
    
    # Function to open file, handling both gzipped and regular files
    def open_file(filename):
        if filename.endswith('.gz'):
            return gzip.open(filename, 'rt')
        else:
            return open(filename, 'r')
    
    print(f"Reading VCF file: {vcf_file}")
    try:
        with open_file(vcf_file) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                    
                fields = line.strip().split('\t')
                if len(fields) < 8:  # Ensure we have at least the INFO field
                    continue
                
                chrom = fields[0]
                pos = int(fields[1])
                rs_id = fields[2]
                ref = fields[3]
                alt = fields[4]
                qual = fields[5]
                filter_status = fields[6]
                info = fields[7]
                
                # Extract required fields from INFO column
                clnsig = extract_info_field(info, "CLNSIG")
                clndn = extract_info_field(info, "CLNDN")
                clnhgvs = extract_info_field(info, "CLNHGVS")
                clnsigconf = extract_info_field(info, "CLNSIGCONF")
                legacy_id = extract_info_field(info, "LEGACY_ID")
                
                # Extract DP from FORMAT and sample columns if available
                dp = ""
                if len(fields) > 9:
                    dp = extract_dp(fields[8], fields[9])
                
                # Extract SnpEff annotations
                ann = extract_info_field(info, "ANN")
                effect, location = extract_snpeff_info(ann)
                
                # Extract LOF information
                lof = extract_info_field(info, "LOF")
                
                # Create variant keys for matching
                # Standard representation and alternative keys for indels
                keys = []
                
                # Handle chromosomes with or without 'chr' prefix
                chrom_variants = [chrom]
                if chrom.startswith('chr'):
                    chrom_variants.append(chrom[3:])
                else:
                    chrom_variants.append(f"chr{chrom}")
                
                # Generate keys for all chromosome variants
                for chr_var in chrom_variants:
                    # Standard key
                    keys.append(f"{chr_var}_{pos}_{ref}_{alt}")
                    
                    # Handle indels
                    if len(ref) > len(alt):  # Deletion
                        for i in range(1, len(ref)):
                            keys.append(f"{chr_var}_{pos+i}_{ref[i]}_{'-'}")
                    elif len(alt) > len(ref):  # Insertion
                        keys.append(f"{chr_var}_{pos}_{'-'}_{alt[len(ref):]}")
                
                variant_info = {
                    "rs": rs_id,
                    "quality": filter_status,
                    "CLNSIG": clnsig,
                    "CLNDN": clndn,
                    "CLNHGVS": clnhgvs,
                    "CLNSIGCONF": clnsigconf,
                    "effect": effect,
                    "location": location,
                    "LOF": lof,
                    "LEGACY_ID": legacy_id,
                    "DP": dp
                }
                
                # Store the variant info for all possible keys
                for key in keys:
                    variants[key] = variant_info
        
        print(f"Parsed {len(variants)} variant positions from VCF")
        return variants
    
    except Exception as e:
        print(f"Error parsing VCF file: {e}", file=sys.stderr)
        sys.exit(1)

def update_tsv(tsv_file, vcf_variants_gatk, vcf_variants_freebayes, output_file, source_column=None):
    """Update the TSV file with information from VCF variants."""
    print(f"Reading TSV file: {tsv_file}")
    try:
        # Read the TSV file
        df = pd.read_csv(tsv_file, sep='\t')
        
        # Try to identify the source column if not specified
        if source_column is None:
            for col in df.columns:
                # Check if column has GATK or FreeBayes values
                col_values = df[col].astype(str)
                if (col_values.str.contains('GATK').any() and 
                    col_values.str.contains('FreeBayes').any()):
                    source_column = col
                    print(f"Found source column: {source_column}")
                    break
            
            # If source column still not found, warn the user
            if source_column is None:
                print("Warning: Could not identify a column containing variant source information (GATK or FreeBayes)")
                print("Will try to match variants from both VCF files.")
        else:
            print(f"Using specified source column: {source_column}")
        
        # Initialize new columns
        new_columns = ["rs", "quality", "CLNSIG", "CLNDN", "CLNHGVS", "CLNSIGCONF", 
                      "effect", "location", "LOF", "LEGACY_ID", "DP"]
        
        for col in new_columns:
            df[col] = ""
        
        # Update rows based on matching variants
        match_count = 0
        gatk_count = 0
        freebayes_count = 0
        
        for idx, row in df.iterrows():
            chrom = str(row['Chr'])
            pos = int(row['Start'])
            ref = str(row['Ref'])
            alt = str(row['Alt'])
            
            # Determine which VCF variants to use based on source
            use_gatk = True
            use_freebayes = True
            
            if source_column:
                source_value = str(row[source_column])
                if "GATK" in source_value:
                    use_gatk = True
                    use_freebayes = False
                elif "FreeBayes" in source_value:
                    use_gatk = False
                    use_freebayes = True
            
            # Try different variant representations for matching
            keys = []
            
            # Handle chromosomes with or without 'chr' prefix
            chrom_variants = [chrom]
            if chrom.startswith('chr'):
                chrom_variants.append(chrom[3:])
            else:
                chrom_variants.append(f"chr{chrom}")
            
            # Generate keys for all chromosome variants
            for chr_var in chrom_variants:
                # Standard key
                keys.append(f"{chr_var}_{pos}_{ref}_{alt}")
                
                # Handle indels
                if alt == '-':  # Deletion in TSV
                    # Try matching with VCF format
                    for prev_base in ['A', 'C', 'G', 'T']:
                        keys.append(f"{chr_var}_{pos-1}_{prev_base+ref}_{prev_base}")
                elif ref == '-':  # Insertion in TSV
                    # Try matching with VCF format
                    for prev_base in ['A', 'C', 'G', 'T']:
                        keys.append(f"{chr_var}_{pos}_{prev_base}_{prev_base+alt}")
            
            # Try all possible keys
            variant_info = None
            
            # Try GATK first if specified
            if use_gatk:
                for key in keys:
                    if key in vcf_variants_gatk:
                        variant_info = vcf_variants_gatk[key]
                        gatk_count += 1
                        break
            
            # If not found in GATK or GATK not specified, try FreeBayes
            if not variant_info and use_freebayes:
                for key in keys:
                    if key in vcf_variants_freebayes:
                        variant_info = vcf_variants_freebayes[key]
                        freebayes_count += 1
                        break
            
            # If found, update the row
            if variant_info:
                for col in new_columns:
                    df.at[idx, col] = variant_info[col]
                match_count += 1
        
        print(f"Matched {match_count} out of {len(df)} variants")
        print(f"GATK matches: {gatk_count}, FreeBayes matches: {freebayes_count}")
        
        # Write the updated dataframe to a new file
        print(f"Writing output to: {output_file}")
        df.to_csv(output_file, sep='\t', index=False)
        print("Done!")
        
    except Exception as e:
        print(f"Error updating TSV file: {e}", file=sys.stderr)
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(description='Add VCF annotations to a TSV file.')
    parser.add_argument('--tsv', required=True, help='Input TSV file')
    parser.add_argument('--gatk_vcf', required=True, help='Input GATK VCF file (gzipped or not)')
    parser.add_argument('--freebayes_vcf', required=True, help='Input FreeBayes VCF file (gzipped or not)')
    parser.add_argument('--output', required=True, help='Output TSV file')
    parser.add_argument('--source_column', help='Name of the column in TSV that contains source information (GATK or FreeBayes)')
    
    args = parser.parse_args()
    
    # Check if input files exist
    if not os.path.exists(args.tsv):
        print(f"Error: TSV file '{args.tsv}' not found", file=sys.stderr)
        sys.exit(1)
    
    if not os.path.exists(args.gatk_vcf):
        print(f"Error: GATK VCF file '{args.gatk_vcf}' not found", file=sys.stderr)
        sys.exit(1)
    
    if not os.path.exists(args.freebayes_vcf):
        print(f"Error: FreeBayes VCF file '{args.freebayes_vcf}' not found", file=sys.stderr)
        sys.exit(1)
    
    # Parse VCF files
    gatk_variants = parse_vcf(args.gatk_vcf)
    freebayes_variants = parse_vcf(args.freebayes_vcf)
    
    # Update TSV file
    update_tsv(args.tsv, gatk_variants, freebayes_variants, args.output, args.source_column)

if __name__ == "__main__":
    main()
