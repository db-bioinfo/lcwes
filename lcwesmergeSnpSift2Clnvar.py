#!/usr/bin/env python3

import csv
import sys

def read_clnvar_file(clnvar_file):
    """Read the clnvar file and create a dictionary of variants with multiple lookup keys."""
    clnvar_variants = {}
    
    with open(clnvar_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            # Store data to add to matched variants
            variant_data = {
                'Filter': row['Filter'],
                'CLNHGVS': row['CLNHGVS'],
                'VAF': row['VAF'],
                'RS': row['RS']  # Added RS field
            }
            
            # Add main keys for lookup
            main_key = (row['Chr'], row['Start'], row['Ref'], row['Alt'], row['Gene'])
            clnvar_variants[main_key] = variant_data
            
            # Add alternate keys for better matching
            # Key without gene for fallback
            no_gene_key = (row['Chr'], row['Start'], row['Ref'], row['Alt'])
            clnvar_variants[no_gene_key] = variant_data
            
            # Add key with End position if different from Start
            if row['End'] != row['Start']:
                end_key = (row['Chr'], row['End'], row['Ref'], row['Alt'], row['Gene'])
                clnvar_variants[end_key] = variant_data
                
                # End position without gene
                end_no_gene_key = (row['Chr'], row['End'], row['Ref'], row['Alt'])
                clnvar_variants[end_no_gene_key] = variant_data
            
            # Handle multiple alternate alleles
            if ',' in row['Alt']:
                alt_alleles = row['Alt'].split(',')
                vaf_values = row['VAF'].split(',') if ',' in row['VAF'] else [row['VAF']] * len(alt_alleles)
                
                for i, alt in enumerate(alt_alleles):
                    if i < len(vaf_values):
                        vaf = vaf_values[i]
                    else:
                        vaf = 'NA'  # In case VAF values don't match alt alleles count
                    
                    # Create individual entries for each alternate allele
                    alt_data = {
                        'Filter': row['Filter'],
                        'CLNHGVS': row['CLNHGVS'],
                        'VAF': vaf,
                        'RS': row['RS']  # Added RS field
                    }
                    
                    alt_key = (row['Chr'], row['Start'], row['Ref'], alt, row['Gene'])
                    clnvar_variants[alt_key] = alt_data
                    
                    alt_no_gene_key = (row['Chr'], row['Start'], row['Ref'], alt)
                    clnvar_variants[alt_no_gene_key] = alt_data
                    
                    if row['End'] != row['Start']:
                        alt_end_key = (row['Chr'], row['End'], row['Ref'], alt, row['Gene'])
                        clnvar_variants[alt_end_key] = alt_data
                        
                        alt_end_no_gene_key = (row['Chr'], row['End'], row['Ref'], alt)
                        clnvar_variants[alt_end_no_gene_key] = alt_data
    
    return clnvar_variants

def process_snpsift_file(snpsift_file, clnvar_variants, output_file):
    """Process the snpsift file and add clnvar information."""
    
    with open(snpsift_file, 'r') as f_in, open(output_file, 'w') as f_out:
        reader = csv.DictReader(f_in, delimiter='\t')
        
        # Create new fieldnames by appending the clnvar columns
        fieldnames = reader.fieldnames + ['Filter', 'CLNHGVS', 'VAF', 'RS']  # Added RS to fieldnames
        
        writer = csv.DictWriter(f_out, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        
        matched_count = 0
        total_count = 0
        
        for row in reader:
            total_count += 1
            # Default values if no match is found
            match_found = False
            filter_val = 'NA'
            clnhgvs_val = 'NA'
            vaf_val = 'NA'
            rs_val = 'NA'  # Default value for RS
            
            # Try multiple matching strategies
            
            # 1. Primary strategy: Try exact match with gene
            key = (row['CHROM'], row['POS'], row['REF'], row['ALT'], row['ANN[0].GENE'])
            if key in clnvar_variants:
                filter_val = clnvar_variants[key]['Filter']
                clnhgvs_val = clnvar_variants[key]['CLNHGVS']
                vaf_val = clnvar_variants[key]['VAF']
                rs_val = clnvar_variants[key]['RS']  # Get RS value
                match_found = True
            
            # 2. Try without gene name
            if not match_found:
                key_no_gene = (row['CHROM'], row['POS'], row['REF'], row['ALT'])
                if key_no_gene in clnvar_variants:
                    filter_val = clnvar_variants[key_no_gene]['Filter']
                    clnhgvs_val = clnvar_variants[key_no_gene]['CLNHGVS']
                    vaf_val = clnvar_variants[key_no_gene]['VAF']
                    rs_val = clnvar_variants[key_no_gene]['RS']  # Get RS value
                    match_found = True
            
            # 3. Try with chromosome standardization if needed
            if not match_found:
                chrom = row['CHROM']
                if not chrom.startswith('chr'):
                    chrom = 'chr' + chrom
                    chr_key = (chrom, row['POS'], row['REF'], row['ALT'], row['ANN[0].GENE'])
                    if chr_key in clnvar_variants:
                        filter_val = clnvar_variants[chr_key]['Filter']
                        clnhgvs_val = clnvar_variants[chr_key]['CLNHGVS']
                        vaf_val = clnvar_variants[chr_key]['VAF']
                        rs_val = clnvar_variants[chr_key]['RS']  # Get RS value
                        match_found = True
                    
                    if not match_found:
                        chr_key_no_gene = (chrom, row['POS'], row['REF'], row['ALT'])
                        if chr_key_no_gene in clnvar_variants:
                            filter_val = clnvar_variants[chr_key_no_gene]['Filter']
                            clnhgvs_val = clnvar_variants[chr_key_no_gene]['CLNHGVS']
                            vaf_val = clnvar_variants[chr_key_no_gene]['VAF']
                            rs_val = clnvar_variants[chr_key_no_gene]['RS']  # Get RS value
                            match_found = True
            
            # 4. For multi-allelic variants, try each allele separately
            if not match_found and ',' in row['ALT']:
                alt_alleles = row['ALT'].split(',')
                for alt in alt_alleles:
                    alt_key = (row['CHROM'], row['POS'], row['REF'], alt, row['ANN[0].GENE'])
                    if alt_key in clnvar_variants:
                        filter_val = clnvar_variants[alt_key]['Filter']
                        clnhgvs_val = clnvar_variants[alt_key]['CLNHGVS']
                        vaf_val = clnvar_variants[alt_key]['VAF']
                        rs_val = clnvar_variants[alt_key]['RS']  # Get RS value
                        match_found = True
                        break
                    
                    alt_key_no_gene = (row['CHROM'], row['POS'], row['REF'], alt)
                    if alt_key_no_gene in clnvar_variants:
                        filter_val = clnvar_variants[alt_key_no_gene]['Filter']
                        clnhgvs_val = clnvar_variants[alt_key_no_gene]['CLNHGVS']
                        vaf_val = clnvar_variants[alt_key_no_gene]['VAF']
                        rs_val = clnvar_variants[alt_key_no_gene]['RS']  # Get RS value
                        match_found = True
                        break
            
            # Add the values to the row and write it out
            row['Filter'] = filter_val
            row['CLNHGVS'] = clnhgvs_val
            row['VAF'] = vaf_val
            row['RS'] = rs_val  # Add RS value to output
            
            if match_found:
                matched_count += 1
            
            writer.writerow(row)
        
        print(f"Total variants processed: {total_count}")
        print(f"Matched variants: {matched_count}")
        print(f"Unmatched variants: {total_count - matched_count}")

def main():
    if len(sys.argv) != 4:
        print("Usage: python script.py <clnvar_file> <snpsift_file> <output_file>")
        sys.exit(1)
    
    clnvar_file = sys.argv[1]
    snpsift_file = sys.argv[2]
    output_file = sys.argv[3]
    
    # Read the clnvar file
    clnvar_variants = read_clnvar_file(clnvar_file)
    
    # Process the snpsift file and add clnvar information
    process_snpsift_file(snpsift_file, clnvar_variants, output_file)
    
    print(f"Merged file saved as {output_file}")

if __name__ == "__main__":
    main()
