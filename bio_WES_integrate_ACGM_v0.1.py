import pandas as pd
import sys
import re

def parse_vcf_header(vcf_content):
    """Extract VCF header and return header lines and the column names."""
    header_lines = []
    for line in vcf_content.split('\n'):
        if line.startswith('#'):
            header_lines.append(line)
            if line.startswith('#CHROM'):
                columns = line[1:].split('\t')
                return header_lines[:-1], columns
    return header_lines, None

def parse_vcf_records(vcf_content):
    """Parse VCF content into a pandas DataFrame."""
    # Skip header lines
    data_lines = [line for line in vcf_content.split('\n') if not line.startswith('#') and line.strip()]
    
    # Parse into DataFrame
    df = pd.DataFrame([line.split('\t') for line in data_lines])
    if len(df.columns) >= 10:
        df.columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'SAMPLE']
    return df

def parse_intervar(intervar_content):
    """Parse InterVar content into a pandas DataFrame."""
    # First, read the header line to get column names
    lines = intervar_content.split('\n')
    header_line = next(line for line in lines if line.startswith('#Chr'))
    
    # Remove the '#' from the header line
    header_line = header_line[1:]
    
    # Split the header into columns
    columns = header_line.split('\t')
    
    # Get data lines (non-header, non-empty lines)
    data_lines = [line.split('\t') for line in lines if not line.startswith('#') and line.strip()]
    
    # Create DataFrame
    df = pd.DataFrame(data_lines, columns=columns)
    
    # Find the InterVar column (it might have spaces in the name)
    intervar_col = [col for col in df.columns if 'InterVar: InterVar and Evidence' in col][0]
    
    # Select and rename relevant columns
    df = df[['Chr', 'Start', 'Ref', 'Alt', intervar_col]]
    df.columns = ['CHROM', 'POS', 'REF', 'ALT', 'INTERVAR']
    
    # Convert POS to string for matching
    df['POS'] = df['POS'].astype(str)
    
    return df

def integrate_intervar_info(vcf_df, intervar_df):
    """Integrate InterVar information into VCF DataFrame."""
    # Create a matching key for both dataframes
    vcf_df['MATCH_KEY'] = vcf_df['CHROM'] + '_' + vcf_df['POS'] + '_' + vcf_df['REF'] + '_' + vcf_df['ALT']
    intervar_df['MATCH_KEY'] = intervar_df['CHROM'] + '_' + intervar_df['POS'] + '_' + intervar_df['REF'] + '_' + intervar_df['ALT']
    
    # Create a dictionary of InterVar info
    intervar_dict = dict(zip(intervar_df['MATCH_KEY'], intervar_df['INTERVAR']))
    
    # Add InterVar info to VCF INFO field
    def add_intervar_info(row):
        match_key = row['MATCH_KEY']
        if match_key in intervar_dict:
            intervar_info = intervar_dict[match_key].strip()
            current_info = row['INFO']
            if not current_info.endswith(';'):
                current_info += ';'
            row['INFO'] = current_info + 'INTERVAR=' + intervar_info
        return row
    
    vcf_df = vcf_df.apply(add_intervar_info, axis=1)
    return vcf_df.drop('MATCH_KEY', axis=1)

def main(vcf_file, intervar_file, output_file):
    """Main function to process and integrate the files."""
    try:
        # Read input files
        print("Reading input files...")
        with open(vcf_file, 'r') as f:
            vcf_content = f.read()
        
        with open(intervar_file, 'r') as f:
            intervar_content = f.read()
        
        # Parse VCF header and content
        print("Parsing VCF file...")
        header_lines, vcf_columns = parse_vcf_header(vcf_content)
        vcf_df = parse_vcf_records(vcf_content)
        
        # Parse InterVar content
        print("Parsing InterVar file...")
        intervar_df = parse_intervar(intervar_content)
        
        # Integrate InterVar information
        print("Integrating InterVar data...")
        result_df = integrate_intervar_info(vcf_df, intervar_df)
        
        # Add INTERVAR to header
        header_lines.append('##INFO=<ID=INTERVAR,Number=.,Type=String,Description="InterVar classification and evidence">')
        
        # Combine header and content
        print("Writing output file...")
        header = '\n'.join(header_lines)
        columns_line = '#' + '\t'.join(vcf_columns)
        content = result_df.to_csv(sep='\t', index=False, header=False)
        
        # Write output
        with open(output_file, 'w') as f:
            f.write(f"{header}\n{columns_line}\n{content}")
            
        print(f"Successfully integrated InterVar data and wrote output to {output_file}")
        
    except Exception as e:
        print(f"Error: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py input.vcf input.intervar.txt output.vcf")
        sys.exit(1)
    
    vcf_file = sys.argv[1]
    intervar_file = sys.argv[2]
    output_file = sys.argv[3]
    
    main(vcf_file, intervar_file, output_file)
