#!/usr/bin/env python3

import pandas as pd
import argparse
import os

def read_variant_file(filename):
    """Read a variant file and return a DataFrame."""
    df = pd.read_csv(filename, sep='\t')
    return df

def create_variant_key(row):
    """Create a unique key for each variant based on Chr, Start, End, Ref, and Alt."""
    return f"{row['Chr']}:{row['Start']}:{row['End']}:{row['Ref']}:{row['Alt']}"

def merge_variants(gatk_file, freebayes_file, output_file):
    """
    Merge variants from GATK and FreeBayes, keeping all GATK variants
    and adding FreeBayes variants that don't exist in GATK.
    """
    print(f"Reading GATK variants from {gatk_file}")
    gatk_df = read_variant_file(gatk_file)
    
    print(f"Reading FreeBayes variants from {freebayes_file}")
    freebayes_df = read_variant_file(freebayes_file)
    
    # Create unique keys for each variant
    gatk_df['variant_key'] = gatk_df.apply(create_variant_key, axis=1)
    freebayes_df['variant_key'] = freebayes_df.apply(create_variant_key, axis=1)
    
    print(f"GATK variants: {len(gatk_df)}")
    print(f"FreeBayes variants: {len(freebayes_df)}")
    
    # Find FreeBayes variants that don't exist in GATK
    gatk_keys = set(gatk_df['variant_key'])
    novel_freebayes = freebayes_df[~freebayes_df['variant_key'].isin(gatk_keys)]
    
    print(f"Unique FreeBayes variants not in GATK: {len(novel_freebayes)}")
    
    # Combine GATK variants with novel FreeBayes variants
    merged_df = pd.concat([gatk_df, novel_freebayes])
    
    # Add a source column to track which caller identified each variant
    merged_df.loc[merged_df['variant_key'].isin(gatk_keys), 'source'] = 'GATK'
    merged_df.loc[~merged_df['variant_key'].isin(gatk_keys), 'source'] = 'FreeBayes'
    
    # Sort by chromosome and position
    merged_df = merged_df.sort_values(['Chr', 'Start'])
    
    # Remove the temporary variant_key column
    merged_df = merged_df.drop('variant_key', axis=1)
    
    # Write the merged variants to file
    print(f"Writing {len(merged_df)} merged variants to {output_file}")
    merged_df.to_csv(output_file, sep='\t', index=False)
    
    print("Summary:")
    print(f"  GATK variants: {len(gatk_df)}")
    print(f"  FreeBayes variants: {len(freebayes_df)}")
    print(f"  Unique FreeBayes variants added: {len(novel_freebayes)}")
    print(f"  Total merged variants: {len(merged_df)}")

def main():
    parser = argparse.ArgumentParser(description='Merge variant calls from GATK and FreeBayes.')
    parser.add_argument('--gatk', required=True, help='GATK variants file')
    parser.add_argument('--freebayes', required=True, help='FreeBayes variants file')
    parser.add_argument('--output', required=True, help='Output merged variants file')
    
    args = parser.parse_args()
    
    if not os.path.exists(args.gatk):
        print(f"Error: GATK file {args.gatk} does not exist")
        return 1
    
    if not os.path.exists(args.freebayes):
        print(f"Error: FreeBayes file {args.freebayes} does not exist")
        return 1
    
    merge_variants(args.gatk, args.freebayes, args.output)
    return 0

if __name__ == "__main__":
    exit(main())
