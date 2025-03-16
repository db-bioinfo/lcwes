#!/usr/bin/env python3
import sys
import csv

def create_variant_key(row):
    """Create a unique key for each variant based on chromosome, position, and alleles"""
    return f"{row[0]}:{row[1]}:{row[2]}:{row[3]}:{row[4]}"

def integrate_variants(gatk_file, freebayes_file, output_file=None):
    """Integrate variants from Freebayes into GATK results"""
    
    try:
        # Read the GATK file
        gatk_variants = []
        gatk_keys = set()
        header = None
        
        with open(gatk_file, 'r') as f:
            reader = csv.reader(f, delimiter='\t')
            header = next(reader)  # Save header
            for row in reader:
                if len(row) >= 5:  # Ensure we have enough columns for the key
                    key = create_variant_key(row)
                    gatk_keys.add(key)
                    gatk_variants.append(row)
                else:
                    print(f"Warning: Skipping malformed GATK row: {row}")
        
        # Read the Freebayes file
        unique_freebayes = []
        chr_counts = {}  # Count variants by chromosome
        
        with open(freebayes_file, 'r') as f:
            reader = csv.reader(f, delimiter='\t')
            fb_header = next(reader)  # Skip header but validate it matches
            
            # Confirm headers match
            if header != fb_header:
                print("Warning: Headers in GATK and Freebayes files don't match exactly.", file=sys.stderr)
                print("Will proceed assuming columns have the same meanings.", file=sys.stderr)
            
            # Get indices for gene and variant type for reporting
            gene_idx = fb_header.index('Ref.Gene') if 'Ref.Gene' in fb_header else None
            type_idx = fb_header.index('Type') if 'Type' in fb_header else None
            func_idx = fb_header.index('Func.refGene') if 'Func.refGene' in fb_header else None
            exonic_func_idx = fb_header.index('ExonicFunc.refGene') if 'ExonicFunc.refGene' in fb_header else None
            
            for row in reader:
                if len(row) >= 5:  # Ensure we have enough columns for the key
                    key = create_variant_key(row)
                    if key not in gatk_keys:
                        unique_freebayes.append(row)
                        
                        # Count variants by chromosome
                        chrom = row[0]
                        chr_counts[chrom] = chr_counts.get(chrom, 0) + 1
                else:
                    print(f"Warning: Skipping malformed Freebayes row: {row}", file=sys.stderr)
        
        # Combine variants
        all_variants = gatk_variants + unique_freebayes
        
        # Sort variants based on final_score (if available)
        final_score_idx = None
        classification_idx = None
        
        if "final_score" in header:
            final_score_idx = header.index("final_score")
        
        if "classification" in header:
            classification_idx = header.index("classification")
            
        if final_score_idx is not None:
            # Sorting by final_score (descending)
            all_variants.sort(key=lambda x: (
                # First sort by classification tier if available
                0 if classification_idx is None or not x[classification_idx].startswith("Tier") 
                else int(x[classification_idx].split()[1]),
                # Then by final score (converting to float, defaulting to 0 if conversion fails)
                -float(x[final_score_idx]) if x[final_score_idx] and x[final_score_idx] != '.' else 0
            ))
        
        # Write the result
        if output_file:
            with open(output_file, 'w') as f:
                writer = csv.writer(f, delimiter='\t')
                writer.writerow(header)
                writer.writerows(all_variants)
        else:
            # Write to stdout
            writer = csv.writer(sys.stdout, delimiter='\t')
            writer.writerow(header)
            writer.writerows(all_variants)
            
        print(f"Original GATK variants: {len(gatk_variants)}", file=sys.stderr)
        print(f"Original Freebayes variants: {len(unique_freebayes) + len(gatk_keys.intersection([create_variant_key(row) for row in unique_freebayes]))}", file=sys.stderr)
        print(f"Unique Freebayes variants found: {len(unique_freebayes)}", file=sys.stderr)
        print(f"Successfully integrated variants. Total variants in output: {len(all_variants)}", file=sys.stderr)
        
        # Print detailed information about integrated variants
        if unique_freebayes:
            print("\nIntegrated variants from Freebayes:", file=sys.stderr)
            print("=================================", file=sys.stderr)
            
            for row in unique_freebayes:
                variant_info = []
                variant_info.append(f"Chr: {row[0]}")
                variant_info.append(f"Pos: {row[1]}-{row[2]}")
                variant_info.append(f"Ref: {row[3]}")
                variant_info.append(f"Alt: {row[4]}")
                
                # Add gene info if available
                if gene_idx is not None and gene_idx < len(row):
                    variant_info.append(f"Gene: {row[gene_idx]}")
                
                # Add variant type info if available
                type_info = []
                if func_idx is not None and func_idx < len(row):
                    type_info.append(row[func_idx])
                if exonic_func_idx is not None and exonic_func_idx < len(row):
                    if row[exonic_func_idx] and row[exonic_func_idx] != '.':
                        type_info.append(row[exonic_func_idx])
                if type_idx is not None and type_idx < len(row):
                    if row[type_idx] and row[type_idx] != '.':
                        type_info.append(row[type_idx])
                
                if type_info:
                    variant_info.append(f"Type: {' - '.join(filter(None, type_info))}")
                
                print(" | ".join(variant_info), file=sys.stderr)
            
            # Print summary by chromosome
            print("\nSummary of integrated variants by chromosome:", file=sys.stderr)
            for chrom in sorted(chr_counts.keys()):
                print(f"{chrom}: {chr_counts[chrom]} variants", file=sys.stderr)
    
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    if len(sys.argv) != 3 and len(sys.argv) != 4:
        print("Usage: python lcmer.py gatk_input freebayes_input [output]")
        sys.exit(1)
    
    gatk_file = sys.argv[1]
    freebayes_file = sys.argv[2]
    output_file = sys.argv[3] if len(sys.argv) == 4 else None
    
    integrate_variants(gatk_file, freebayes_file, output_file)
