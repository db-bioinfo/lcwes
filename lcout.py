import csv
import argparse
import sys
import os

def integrate_variants(gatk_file, freebayes_file, output_file):
    """
    Simple integration of Freebayes variants into GATK by line number.
    If variant X is at line 5 in Freebayes and doesn't exist in GATK, 
    insert it at line 5 in the output and shift all subsequent GATK variants down.
    
    Args:
        gatk_file: Path to GATK variants TSV file
        freebayes_file: Path to Freebayes variants TSV file
        output_file: Path to write the integrated variants
    """
    try:
        # Read the GATK file
        gatk_variants = []
        gatk_ids = set()
        
        with open(gatk_file, 'r') as f:
            reader = csv.reader(f, delimiter='\t')
            header = next(reader)
            for row in reader:
                gatk_variants.append(row)
                # Create a unique identifier using Chr_Start_Ref_Alt
                variant_id = f"{row[0]}_{row[1]}_{row[2]}_{row[3]}"
                gatk_ids.add(variant_id)
        
        print(f"Read {len(gatk_variants)} variants from GATK")
        
        # Read the Freebayes file
        freebayes_variants = []
        
        with open(freebayes_file, 'r') as f:
            reader = csv.reader(f, delimiter='\t')
            freebayes_header = next(reader)  # Skip header
            
            for row in reader:
                # Create a unique identifier using Chr_Start_Ref_Alt
                variant_id = f"{row[0]}_{row[1]}_{row[2]}_{row[3]}"
                # Add the variant to our list with its ID for quick lookup
                freebayes_variants.append((variant_id, row))
        
        print(f"Read {len(freebayes_variants)} variants from Freebayes")
        
        # Create the result array, starting with the GATK variants
        result_variants = list(gatk_variants)
        
        # For each Freebayes variant, check if it's already in GATK
        # If not, insert it at the same position in the result
        inserted_count = 0
        
        for idx, (variant_id, row) in enumerate(freebayes_variants):
            if variant_id not in gatk_ids:
                # Insert at the same line number (adjusted for already inserted variants)
                insert_pos = min(idx + inserted_count, len(result_variants))
                result_variants.insert(insert_pos, row)
                inserted_count += 1
                
                # Log the insertion
                variant_name = f"{row[4]}:{row[0]}:{row[1]}" if len(row) > 4 else f"{row[0]}:{row[1]}"
                print(f"Inserted {variant_name} at position {insert_pos+2} (from line {idx+2} in Freebayes)")
        
        print(f"Inserted {inserted_count} unique variants from Freebayes")
        
        # Write to output file
        with open(output_file, 'w', newline='') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow(header)
            writer.writerows(result_variants)
        
        print(f"Wrote {len(result_variants)} variants to {output_file}")
        
        return True
    except Exception as e:
        print(f"Error integrating variants: {e}")
        import traceback
        traceback.print_exc()
        return False

def main():
    parser = argparse.ArgumentParser(description='Integrate Freebayes variants into GATK output')
    parser.add_argument('--gatk', required=True, help='Path to GATK variants TSV file')
    parser.add_argument('--freebayes', required=True, help='Path to Freebayes variants TSV file')
    parser.add_argument('--output', required=True, help='Path to write integrated variants')
    
    args = parser.parse_args()
    
    # Check if input files exist
    if not os.path.exists(args.gatk):
        print(f"GATK file not found: {args.gatk}")
        sys.exit(1)
    
    if not os.path.exists(args.freebayes):
        print(f"Freebayes file not found: {args.freebayes}")
        sys.exit(1)
    
    # Create output directory if it doesn't exist
    output_dir = os.path.dirname(args.output)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Integrate the variants
    success = integrate_variants(args.gatk, args.freebayes, args.output)
    
    if not success:
        sys.exit(1)

if __name__ == "__main__":
    main()
