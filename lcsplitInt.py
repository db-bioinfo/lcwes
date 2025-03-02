#!/usr/bin/env python3
import re
import sys
import csv

def extract_classification(intervar_str):
    """Extract the ACMG classification from the InterVar string."""
    # Try to match text between "InterVar: " and the first evidence criterion
    classification_match = re.search(r'InterVar: ([A-Za-z ]+?)(?=PVS1|$)', intervar_str)
    if classification_match and classification_match.group(1).strip():
        return classification_match.group(1).strip()
    
    # If that fails, try to get any text after "InterVar: "
    alt_match = re.search(r'InterVar: ([A-Za-z ]+)', intervar_str)
    return alt_match.group(1).strip() if alt_match else ""

def extract_met_criteria(intervar_str):
    """Extract the met ACMG criteria from the InterVar string."""
    met_criteria = []
    
    # Extract PVS1
    if 'PVS1=1' in intervar_str:
        met_criteria.append('PVS1')
    
    # Extract PS array
    ps_match = re.search(r'PS=\[([\d, ]+)\]', intervar_str)
    if ps_match:
        ps_array = [s.strip() for s in ps_match.group(1).split(',')]
        for i, value in enumerate(ps_array):
            if value == '1':
                met_criteria.append(f'PS{i+1}')
    
    # Extract PM array
    pm_match = re.search(r'PM=\[([\d, ]+)\]', intervar_str)
    if pm_match:
        pm_array = [s.strip() for s in pm_match.group(1).split(',')]
        for i, value in enumerate(pm_array):
            if value == '1':
                met_criteria.append(f'PM{i+1}')
    
    # Extract PP array
    pp_match = re.search(r'PP=\[([\d, ]+)\]', intervar_str)
    if pp_match:
        pp_array = [s.strip() for s in pp_match.group(1).split(',')]
        for i, value in enumerate(pp_array):
            if value == '1':
                met_criteria.append(f'PP{i+1}')
    
    return ', '.join(met_criteria)

def process_file(input_file, output_file):
    """Process the TSV file to split the InterVar column."""
    with open(input_file, 'r', newline='') as infile, open(output_file, 'w', newline='') as outfile:
        reader = csv.reader(infile, delimiter='\t')
        writer = csv.writer(outfile, delimiter='\t')
        
        # Process header
        headers = next(reader)
        intervar_column_index = -1
        
        for i, header in enumerate(headers):
            if header.strip() == 'InterVar: InterVar and Evidence':
                intervar_column_index = i
                break
        
        if intervar_column_index == -1:
            print("Error: 'InterVar: InterVar and Evidence' column not found in the file.")
            return
        
        # Create new headers
        new_headers = (
            headers[:intervar_column_index] + 
            ['ACMG_Classification', 'ACMG_Met_Criteria'] + 
            headers[intervar_column_index+1:]
        )
        writer.writerow(new_headers)
        
        # Process data rows
        for row in reader:
            if not row or len(row) <= intervar_column_index:
                continue  # Skip empty or malformed rows
            
            intervar_str = row[intervar_column_index]
            classification = extract_classification(intervar_str)
            met_criteria = extract_met_criteria(intervar_str)
            
            new_row = (
                row[:intervar_column_index] + 
                [classification, met_criteria] + 
                row[intervar_column_index+1:]
            )
            writer.writerow(new_row)

def main():
    if len(sys.argv) != 3:
        print("Usage: python script.py input.tsv output.tsv")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    process_file(input_file, output_file)
    print(f"Processed {input_file} and saved to {output_file}")

if __name__ == "__main__":
    main()
