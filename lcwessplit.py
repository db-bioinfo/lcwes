#!/usr/bin/env python3

import pandas as pd
import re
import sys

def replace_intervar_column(input_file, output_file):
    """
    Replace the 'InterVar: InterVar and Evidence' column with two new columns:
    'ACMG' and 'ACMG Criteria'
    
    Args:
        input_file (str): Path to the input TSV file
        output_file (str): Path to the output TSV file
    """
    print(f"Reading {input_file}...")
    
    # Read the TSV file
    try:
        df = pd.read_csv(input_file, sep='\t')
        print(f"Successfully read file with {len(df)} rows and {len(df.columns)} columns")
    except Exception as e:
        print(f"Error reading input file: {e}")
        sys.exit(1)
    
    # Check if the InterVar column exists
    intervar_col = ' InterVar: InterVar and Evidence '
    if intervar_col not in df.columns:
        print(f"Column '{intervar_col}' not found in the input file.")
        print("Available columns:")
        for col in df.columns:
            print(f"  - '{col}'")
        sys.exit(1)
    
    print(f"Processing the InterVar column: '{intervar_col}'")
    
    # Get the position of the InterVar column
    col_pos = list(df.columns).index(intervar_col)
    
    # Create temporary columns to store the extracted values
    df['ACMG_temp'] = ""
    df['ACMG Criteria_temp'] = ""
    
    # Process each row
    for idx, row in df.iterrows():
        intervar_value = str(row[intervar_col])
        
        # Skip if the value is empty or NaN
        if not intervar_value or intervar_value == 'nan':
            continue
        
        # Extract ACMG classification
        acmg_match = re.search(r'InterVar: (.*?) PVS', intervar_value)
        if acmg_match:
            df.at[idx, 'ACMG_temp'] = acmg_match.group(1).strip()
        
        # Extract criteria
        criteria = []
        
        # PVS1
        pvs_match = re.search(r'PVS1=(\d+)', intervar_value)
        if pvs_match and pvs_match.group(1) != '0':
            criteria.append("PVS1")
        
        # PS criteria
        ps_match = re.search(r'PS=\[(.*?)\]', intervar_value)
        if ps_match:
            ps_values = ps_match.group(1).split(', ')
            for i, val in enumerate(ps_values, 1):
                if val != '0':
                    criteria.append(f"PS{i}")
        
        # PM criteria
        pm_match = re.search(r'PM=\[(.*?)\]', intervar_value)
        if pm_match:
            pm_values = pm_match.group(1).split(', ')
            for i, val in enumerate(pm_values, 1):
                if val != '0':
                    criteria.append(f"PM{i}")
        
        # PP criteria
        pp_match = re.search(r'PP=\[(.*?)\]', intervar_value)
        if pp_match:
            pp_values = pp_match.group(1).split(', ')
            for i, val in enumerate(pp_values, 1):
                if val != '0':
                    criteria.append(f"PP{i}")
        
        df.at[idx, 'ACMG Criteria_temp'] = ", ".join(criteria)
    
    # Remove the original InterVar column
    df = df.drop(columns=[intervar_col])
    
    # Get list of all columns
    cols = list(df.columns)
    
    # Remove temporary columns from the list
    cols.remove('ACMG_temp')
    cols.remove('ACMG Criteria_temp')
    
    # Create a new list with columns in the right order
    new_cols = cols[:col_pos] + ['ACMG', 'ACMG Criteria'] + cols[col_pos:]
    
    # Create a new DataFrame with columns in the right order
    new_df = pd.DataFrame()
    
    # Add all columns except the temporary ones
    for col in cols:
        new_df[col] = df[col]
    
    # Add the new columns at the right position
    new_df.insert(col_pos, 'ACMG', df['ACMG_temp'])
    new_df.insert(col_pos + 1, 'ACMG Criteria', df['ACMG Criteria_temp'])
    
    print(f"Replaced column '{intervar_col}' with 'ACMG' and 'ACMG Criteria'")
    
    # Write output
    try:
        new_df.to_csv(output_file, sep='\t', index=False)
        print(f"Successfully wrote output to {output_file}")
    except Exception as e:
        print(f"Error writing output file: {e}")
        sys.exit(1)

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python replace_intervar.py input.tsv output.tsv")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    replace_intervar_column(input_file, output_file)
