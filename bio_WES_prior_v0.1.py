import pandas as pd
from cyvcf2 import VCF
from collections import defaultdict
import sys
import os

class ClinicalVariantPrioritizer:
    def __init__(self):
        # ACMG evidence weights
        self.acmg_evidence_scores = {
            'PVS1': 350,  # Very Strong evidence of pathogenicity
            'PS': 250,    # Strong evidence of pathogenicity
            'PM': 150,    # Moderate evidence of pathogenicity
            'PP': 50,     # Supporting evidence of pathogenicity
            'BA1': -350,  # Stand-alone evidence of benign impact
            'BS': -250,   # Strong evidence of benign impact
            'BP': -50     # Supporting evidence of benign impact
        }
        
        # ACMG classification weights (additional to evidence-based score)
        self.acmg_classification_scores = {
            'Pathogenic': 500,
            'Likely pathogenic': 400,
            'Uncertain significance': 0,
            'Likely benign': -400,
            'Benign': -500
        }
        
        # Clinical significance categories and their weights
        self.clinical_significance_scores = {
            'pathogenic': 1000,
            'likely_pathogenic': 800,
            'uncertain_significance': 400,
            'likely_benign': -500,
            'benign': -800,
            'not_provided': 0,
            'conflicting': 200
        }
        
        # Variant impact categories and their weights
        self.impact_scores = {
            'HIGH': 500,
            'MODERATE': 300,
            'LOW': 100,
            'MODIFIER': 0
        }
        
        # Specific variant consequences and their clinical relevance scores
        self.consequence_scores = {
            'frameshift': 500,
            'stop_gained': 500,
            'stop_lost': 450,
            'start_lost': 450,
            'splice_acceptor': 450,
            'splice_donor': 450,
            'missense': 300,
            'splice_region': 200,
            'inframe': 200,
            'synonymous': 50,
            'intron': 10,
            '5_prime_utr': 20,
            '3_prime_utr': 20,
            'upstream': 0,
            'downstream': 0
        }

    def parse_intervar_info(self, info_str):
        """Parse InterVar information from INFO field."""
        if not info_str or 'InterVar:' not in info_str:
            return "Not available", "Not available"
            
        # Extract classification
        classification = "Not available"
        criteria_met = "Not available"
        
        try:
            # Get the classification
            if "InterVar:" in info_str:
                classification = info_str.split("InterVar:")[1].split("PVS")[0].strip()
            
            # Extract criteria details
            criteria_parts = []
            
            # Extract PVS
            if "PVS1=" in info_str:
                pvs_val = info_str.split("PVS1=")[1].split()[0]
                if pvs_val != "0":
                    criteria_parts.append(f"PVS1")
            
            # Extract PS array
            if "PS=[" in info_str:
                ps_arr = info_str.split("PS=[")[1].split("]")[0].split(",")
                for i, val in enumerate(ps_arr):
                    if val.strip() != "0":
                        criteria_parts.append(f"PS{i+1}")
            
            # Extract PM array
            if "PM=[" in info_str:
                pm_arr = info_str.split("PM=[")[1].split("]")[0].split(",")
                for i, val in enumerate(pm_arr):
                    if val.strip() != "0":
                        criteria_parts.append(f"PM{i+1}")
            
            # Extract PP array
            if "PP=[" in info_str:
                pp_arr = info_str.split("PP=[")[1].split("]")[0].split(",")
                for i, val in enumerate(pp_arr):
                    if val.strip() != "0":
                        criteria_parts.append(f"PP{i+1}")
            
            # Extract BA1
            if "BA1=" in info_str:
                ba1_val = info_str.split("BA1=")[1].split()[0]
                if ba1_val != "0":
                    criteria_parts.append("BA1")
            
            # Extract BS array
            if "BS=[" in info_str:
                bs_arr = info_str.split("BS=[")[1].split("]")[0].split(",")
                for i, val in enumerate(bs_arr):
                    if val.strip() != "0":
                        criteria_parts.append(f"BS{i+1}")
            
            # Extract BP array
            if "BP=[" in info_str:
                bp_arr = info_str.split("BP=[")[1].split("]")[0].split(",")
                for i, val in enumerate(bp_arr):
                    if val.strip() != "0":
                        criteria_parts.append(f"BP{i+1}")
            
            # Join all met criteria
            criteria_met = ", ".join(criteria_parts) if criteria_parts else "None"
            
        except Exception as e:
            print(f"Error parsing InterVar info: {e}")
            
        return classification, criteria_met

    def get_location_info(self, effect_type, location_field, hgvs_c):
        """Extract standardized location information."""
        location = "Unknown"
        effect_lower = effect_type.lower()

        try:
            # Coding variants that are always in exons
            exonic_effects = {
                'frameshift', 'stop_gained', 'stop_lost', 'start_lost',
                'missense', 'synonymous', 'coding', 'inframe',
                'protein_altering', 'initiator_codon', 'incomplete_terminal_codon'
            }
            
            # Check if it's a coding variant
            is_coding_variant = any(effect in effect_lower for effect in exonic_effects)
            
            if is_coding_variant:
                if location_field and '/' in location_field:
                    current, total = location_field.split('/')
                    location = f"Exon {current}/{total}"
                else:
                    location = "Exonic"
                    
            elif 'intron' in effect_lower:
                if location_field and '/' in location_field:
                    current, total = location_field.split('/')
                    location = f"Intron {current}/{total}"
                else:
                    location = "Intronic"
                    
            elif 'splice' in effect_lower:
                if 'acceptor' in effect_lower:
                    location = "Splice Acceptor"
                elif 'donor' in effect_lower:
                    location = "Splice Donor"
                else:
                    location = "Splice Region"
                    
            elif '5_prime_utr' in effect_lower:
                location = "5' UTR"
            elif '3_prime_utr' in effect_lower:
                location = "3' UTR"
            elif 'upstream' in effect_lower:
                location = "Upstream"
            elif 'downstream' in effect_lower:
                location = "Downstream"
            elif 'noncoding' in effect_lower or 'non_coding' in effect_lower:
                if location_field and '/' in location_field:
                    current, total = location_field.split('/')
                    location = f"Non-coding exon {current}/{total}"
                else:
                    location = "Non-coding"
            else:
                location = "Other"
                
        except Exception as e:
            location = "Unknown"

        return location

    def parse_ann_field(self, ann_string):
        """Parse the ANN field from SnpEff annotation."""
        effects = ann_string.split(',')
        all_effects = []
        for effect in effects:
            fields = effect.split('|')
            if len(fields) >= 10:
                # Get location information using the new method
                location = self.get_location_info(
                    fields[1],  # effect type
                    fields[8] if len(fields) > 8 else '',  # location field
                    fields[9] if len(fields) > 9 else ''   # HGVS coding
                )
                
                effect_dict = {
                    'allele': fields[0],
                    'effect': fields[1],
                    'impact': fields[2],
                    'gene': fields[3],
                    'gene_id': fields[4],
                    'feature_type': fields[5],
                    'feature_id': fields[6],
                    'transcript_biotype': fields[7],
                    'location': location,
                    'hgvs_c': fields[9] if fields[9] != '' else None,
                    'hgvs_p': fields[10] if len(fields) > 10 and fields[10] != '' else None
                }
                all_effects.append(effect_dict)
        
        # Get the most severe effect based on impact hierarchy
        most_severe = max(all_effects, key=lambda x: self.impact_scores.get(x['impact'], 0))
        
        # Format transcript with HGVS notation if available
        if most_severe['feature_id']:
            notation_parts = []
            notation_parts.append(most_severe['feature_id'])
            if most_severe['hgvs_c']:
                notation_parts.append(most_severe['hgvs_c'])
            if most_severe['hgvs_p']:
                notation_parts.append(most_severe['hgvs_p'])
            
            most_severe['feature_id'] = ':'.join(filter(None, notation_parts))
        
        return most_severe

    def get_clinical_score(self, clinsig_str):
        """Calculate clinical significance score."""
        if not clinsig_str or clinsig_str == 'Not available':
            return 0
            
        score = 0
        clinsig_str = str(clinsig_str).lower()
        
        # Handle multiple clinical significance annotations
        for significance, value in self.clinical_significance_scores.items():
            if significance in clinsig_str:
                score = max(score, value)
        
        return score

    def get_consequence_score(self, effect):
        """Calculate consequence score based on variant effect."""
        score = 0
        effect_lower = str(effect).lower()
        
        for consequence, value in self.consequence_scores.items():
            if consequence in effect_lower:
                score = max(score, value)
                
        return score

    def calculate_acmg_score(self, intervar_info):
        """Calculate score based on ACMG criteria and classification."""
        score = 0
        
        if not intervar_info or 'InterVar:' not in intervar_info:
            return score
            
        # Extract classification and add classification-based score
        classification = intervar_info.split("InterVar:")[1].split()[0]
        if classification in self.acmg_classification_scores:
            score += self.acmg_classification_scores[classification]
        
        # Add scores for individual criteria
        # PVS1
        if "PVS1=1" in intervar_info:
            score += self.acmg_evidence_scores['PVS1']
            
        # PS array
        if "PS=[" in intervar_info:
            ps_arr = intervar_info.split("PS=[")[1].split("]")[0].split(",")
            score += sum(1 for x in ps_arr if x.strip() == "1") * self.acmg_evidence_scores['PS']
            
        # PM array
        if "PM=[" in intervar_info:
            pm_arr = intervar_info.split("PM=[")[1].split("]")[0].split(",")
            score += sum(1 for x in pm_arr if x.strip() == "1") * self.acmg_evidence_scores['PM']
            
        # PP array
        if "PP=[" in intervar_info:
            pp_arr = intervar_info.split("PP=[")[1].split("]")[0].split(",")
            score += sum(1 for x in pp_arr if x.strip() == "1") * self.acmg_evidence_scores['PP']
            
        # BA1
        if "BA1=1" in intervar_info:
            score += self.acmg_evidence_scores['BA1']
            
        # BS array
        if "BS=[" in intervar_info:
            bs_arr = intervar_info.split("BS=[")[1].split("]")[0].split(",")
            score += sum(1 for x in bs_arr if x.strip() == "1") * self.acmg_evidence_scores['BS']
            
        # BP array
        if "BP=[" in intervar_info:
            bp_arr = intervar_info.split("BP=[")[1].split("]")[0].split(",")
            score += sum(1 for x in bp_arr if x.strip() == "1") * self.acmg_evidence_scores['BP']
            
        return score

    def calculate_quality_score(self, qual, depth):
        """Calculate quality score based on QUAL and DP fields."""
        quality_score = 0
        
        # Quality score (max 100)
        if qual is not None:
            quality_score += min(float(qual), 1000) / 10
            
        # Depth score (max 100)
        if depth:
            try:
                dp = int(depth)
                if dp >= 20:  # Minimum clinical coverage
                    quality_score += min(dp, 100)
                else:
                    quality_score -= (20 - dp) * 10  # Penalize low coverage
            except (ValueError, TypeError):
                pass
                
        return quality_score

    def prioritize_variants(self, vcf_file, output_file, max_rows=25000):
        """Main function to prioritize variants and export to Excel."""
        variants = []
        vcf_reader = VCF(vcf_file)
        
        for record in vcf_reader:
            # Parse annotations
            ann_dict = self.parse_ann_field(record.INFO.get('ANN', ''))
            
            # Get clinical significance
            clinsig = record.INFO.get('CLNSIG', 'Not available')
            
            # Get disease name and review status
            clndn = record.INFO.get('CLNDN', 'Not available')
            clnrevstat = record.INFO.get('CLNREVSTAT', 'Not available')
            
            # Get allele frequencies
            af = record.INFO.get('AF', 'Not available')
            gnomad_af = record.INFO.get('gnomAD_AF', 'Not available')
            
            # Get FILTER field
            filter_value = 'PASS' if record.FILTER is None else ','.join(record.FILTER)
            
            # Get COSMIC Legacy ID
            cosmic_id = record.INFO.get('LEGACY_ID', 'Not available')
            
            # Parse InterVar information
            acmg_class, acmg_criteria = self.parse_intervar_info(record.INFO.get('INTERVAR', ''))
            
            # Calculate scores
            clinical_score = self.get_clinical_score(clinsig)
            consequence_score = self.get_consequence_score(ann_dict['effect'])
            impact_score = self.impact_scores.get(ann_dict['impact'], 0)
            quality_score = self.calculate_quality_score(record.QUAL, record.INFO.get('DP', None))
            
            # Calculate ACMG score
            acmg_score = self.calculate_acmg_score(record.INFO.get('INTERVAR', ''))
            
            # Total priority score (now including ACMG score)
            priority_score = clinical_score + consequence_score + impact_score + quality_score + acmg_score
            
            # Generate ClinVar link if needed
            clinvar_link = ''
            clinsig_lower = str(clinsig).lower()
            if ('pathogenic' in clinsig_lower or 
                'likely_pathogenic' in clinsig_lower or 
                'conflicting' in clinsig_lower):
                allele_id = record.INFO.get('ALLELEID', '')
                if allele_id:
                    clinvar_link = f"https://www.ncbi.nlm.nih.gov/clinvar/variation/{allele_id}/"
            
            variant = {
                'CHROM': record.CHROM,
                'POS': record.POS,
                'REF': record.REF,
                'ALT': ','.join(str(alt) for alt in record.ALT),
                'Gene': ann_dict['gene'],
                'dbSNP_ID': record.ID if record.ID else 'Novel',
                'COSMIC_ID': cosmic_id,
                'FILTER': filter_value,
                'Clinical_Significance': clinsig,
                'ClinVar_Link': clinvar_link,
                'CLNHGVS': record.INFO.get('CLNHGVS', 'Not available'),
                'Disease': clndn,
                'AF': af,
                'gnomAD_AF': gnomad_af,
                'Effect': ann_dict['effect'],
                'Impact': ann_dict['impact'],
                'Location': ann_dict['location'],
                'Transcript': ann_dict['feature_id'],
                'Quality': record.QUAL,
                'Depth': record.INFO.get('DP', ''),
                'Review_Status': clnrevstat,
                'ACMG_Classification': acmg_class,
                'ACMG_Met_Criteria': acmg_criteria,
                'Priority_Score': priority_score
            }
            
            variants.append(variant)
        
        # Convert to DataFrame with specified column order
        df = pd.DataFrame(variants)
        column_order = [
            'CHROM', 'POS', 'REF', 'ALT', 'Gene', 'Location', 'dbSNP_ID', 'COSMIC_ID', 'Transcript', 'FILTER',
            'Clinical_Significance', 'ACMG_Classification', 'ACMG_Met_Criteria',
            'ClinVar_Link', 'CLNHGVS', 'Disease', 'AF', 'gnomAD_AF', 'Effect', 'Impact', 'Quality', 'Depth',
            'Review_Status', 'Priority_Score'
        ]
        df_sorted = df[column_order].sort_values('Priority_Score', ascending=False)
        
        # Limit the number of rows to max_rows
        df_sorted = df_sorted.head(max_rows)
        
        # Export to Excel with formatting
        writer = pd.ExcelWriter(output_file, engine='xlsxwriter')
        df_sorted.to_excel(writer, sheet_name='Prioritized Variants', index=False)
        
        # Get workbook and worksheet objects
        workbook = writer.book
        worksheet = writer.sheets['Prioritized Variants']
        
        # Add formats
        header_format = workbook.add_format({
            'bold': True,
            'text_wrap': True,
            'valign': 'top',
            'bg_color': '#D9E1F2',
            'border': 1
        })
        
        # Color formats for ACMG classifications
        acmg_path_format = workbook.add_format({'bg_color': '#FF0000'})      # Red for pathogenic
        acmg_lpath_format = workbook.add_format({'bg_color': '#FF9999'})     # Light red for likely pathogenic
        acmg_vus_format = workbook.add_format({'bg_color': '#FFEB9C'})       # Yellow for VUS
        acmg_lben_format = workbook.add_format({'bg_color': '#C6EFCE'})      # Light green for likely benign
        acmg_ben_format = workbook.add_format({'bg_color': '#92D050'})       # Green for benign
        
        # Updated color formats for ClinVar
        path_format = workbook.add_format({'bg_color': '#FF0000'})     # Bright red for pathogenic
        lpath_format = workbook.add_format({'bg_color': '#FF0000'})    # Same red for likely pathogenic
        vus_format = workbook.add_format({'bg_color': '#FFEB9C'})      # Yellow for VUS
        benign_format = workbook.add_format({'bg_color': '#C6EFCE'})   # Green for benign
        conflict_format = workbook.add_format({'bg_color': '#FFA07A'}) # Light salmon for conflicting
        
        # FILTER formats
        filter_pass_format = workbook.add_format({'bg_color': '#E2EFDA'})  # Light green for PASS
        filter_fail_format = workbook.add_format({'bg_color': '#FFC7CE'})  # Light red for non-PASS
        
        # Link format
        link_format = workbook.add_format({
            'color': 'blue',
            'underline': True
        })
        
        # Apply column formats
        for col_num, value in enumerate(df_sorted.columns.values):
            worksheet.write(0, col_num, value, header_format)
            # Make columns wider for better readability
            if value in ['Disease', 'Review_Status', 'Clinical_Significance', 'FILTER', 'Location', 
                        'ClinVar_Link', 'CLNHGVS', 'ACMG_Classification', 'ACMG_Met_Criteria', 'COSMIC_ID']:
                worksheet.set_column(col_num, col_num, 30)
            else:
                worksheet.set_column(col_num, col_num, 15)
        
        # Apply conditional formatting for clinical significance, ACMG, and write links
        clin_sig_col = df_sorted.columns.get_loc('Clinical_Significance')
        acmg_class_col = df_sorted.columns.get_loc('ACMG_Classification')
        filter_col = df_sorted.columns.get_loc('FILTER')
        link_col = df_sorted.columns.get_loc('ClinVar_Link')
        
        for row_num, (clin_sig, acmg_class, filter_val, link) in enumerate(
            zip(df_sorted['Clinical_Significance'], 
                df_sorted['ACMG_Classification'],
                df_sorted['FILTER'], 
                df_sorted['ClinVar_Link']), start=1):
            
            # Format Clinical Significance column
            clin_sig_lower = str(clin_sig).lower()
            if 'conflicting' in clin_sig_lower:
                worksheet.write(row_num, clin_sig_col, clin_sig, conflict_format)
            elif 'likely_pathogenic' in clin_sig_lower:
                worksheet.write(row_num, clin_sig_col, clin_sig, lpath_format)
            elif 'pathogenic' in clin_sig_lower:
                worksheet.write(row_num, clin_sig_col, clin_sig, path_format)
            elif 'uncertain_significance' in clin_sig_lower:
                worksheet.write(row_num, clin_sig_col, clin_sig, vus_format)
            elif 'benign' in clin_sig_lower:
                worksheet.write(row_num, clin_sig_col, clin_sig, benign_format)
            
            # Format ACMG Classification column
            acmg_lower = str(acmg_class).lower()
            if 'pathogenic' in acmg_lower and 'likely' not in acmg_lower:
                worksheet.write(row_num, acmg_class_col, acmg_class, acmg_path_format)
            elif 'likely pathogenic' in acmg_lower:
                worksheet.write(row_num, acmg_class_col, acmg_class, acmg_lpath_format)
            elif 'uncertain significance' in acmg_lower:
                worksheet.write(row_num, acmg_class_col, acmg_class, acmg_vus_format)
            elif 'likely benign' in acmg_lower:
                worksheet.write(row_num, acmg_class_col, acmg_class, acmg_lben_format)
            elif 'benign' in acmg_lower:
                worksheet.write(row_num, acmg_class_col, acmg_class, acmg_ben_format)
                
            # Format FILTER column
            if filter_val == 'PASS':
                worksheet.write(row_num, filter_col, filter_val, filter_pass_format)
            else:
                worksheet.write(row_num, filter_col, filter_val, filter_fail_format)
            
            # Write ClinVar link as hyperlink if available
            if link:
                worksheet.write_url(row_num, link_col, link, link_format, link)
            else:
                worksheet.write(row_num, link_col, '')
        
        writer.close()
        return df_sorted

def main():
    # Check if correct number of arguments provided
    if len(sys.argv) != 3:
        print("Usage: python script.py input.vcf output.xlsx")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    # Check if input file exists
    if not os.path.exists(input_file):
        print(f"Error: Input file '{input_file}' does not exist")
        sys.exit(1)
    
    # Check if input file is a VCF
    if not input_file.lower().endswith('.vcf'):
        print("Warning: Input file does not have .vcf extension")
    
    # Check if output file has xlsx extension
    if not output_file.lower().endswith('.xlsx'):
        output_file += '.xlsx'
    
    try:
        prioritizer = ClinicalVariantPrioritizer()
        prioritized_df = prioritizer.prioritize_variants(input_file, output_file)
        print(f"Successfully prioritized variants. Results saved to {output_file}")
        if len(prioritized_df) >= 1048576:
            print("Note: Output has been limited to the top 1,048,576 variants by priority score")
    except Exception as e:
        print(f"Error occurred: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main()
