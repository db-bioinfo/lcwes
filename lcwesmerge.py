#!/usr/bin/env python3

import pandas as pd
import sys
import logging
from typing import List, Tuple, Dict
from datetime import datetime

# Set up logging with timestamp in filename
log_filename = f"variant_merge_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(log_filename),
        logging.StreamHandler()
    ]
)

class VariantMerger:
    """
    Enhanced variant merger for clinical annotation data.
    Handles InterVar and MultiAnno file formats with detailed validation.
    """
    
    def __init__(self, intervar_file: str, multianno_file: str):
        self.intervar_file = intervar_file
        self.multianno_file = multianno_file
        self.multianno_columns = [
            'Polyphen2_HDIV_pred', 'MutationTaster_pred', 'REVEL_score',
            'ClinPred_score', 'AlphaMissense_pred', 'Aloft_pred',
            'PrimateAI_pred', 'BayesDel_addAF_pred', 'CLNDN',
            'phyloP100way_vertebrate', 'GERP++_RS'
        ]
        self.stats = {}

    def create_variant_key(self, df: pd.DataFrame) -> pd.Series:
        """
        Creates a unique key for variant matching.
        Format: chr_start_end_ref_alt
        """
        return df['Chr'].astype(str) + '_' + \
               df['Start'].astype(str) + '_' + \
               df['End'].astype(str) + '_' + \
               df['Ref'].astype(str) + '_' + \
               df['Alt'].astype(str)

    def validate_dataframe(self, df: pd.DataFrame, name: str) -> None:
        """
        Performs comprehensive validation of input dataframe.
        Analyzes duplicates, missing values, and gene patterns.
        """
        # Check for missing values in key columns
        key_cols = ['Chr', 'Start', 'End', 'Ref', 'Alt']
        missing_key_vals = df[key_cols].isnull().sum()
        
        if missing_key_vals.sum() > 0:
            for col in key_cols:
                if missing_key_vals[col] > 0:
                    logging.warning(f"{name}: Found {missing_key_vals[col]} missing values in {col}")
        
        # Enhanced duplicate analysis
        variant_keys = self.create_variant_key(df)
        duplicates = variant_keys.duplicated(keep='first')
        num_duplicates = duplicates.sum()
        
        if num_duplicates > 0:
            # Analyze duplicate patterns
            duplicate_variants = df[variant_keys.isin(variant_keys[duplicates])]
            duplicate_counts = variant_keys[variant_keys.isin(variant_keys[duplicates])].value_counts()
            
            # Store detailed duplicate statistics
            self.stats[f'{name}_duplicates'] = num_duplicates
            self.stats[f'{name}_unique_duplicated_variants'] = len(duplicate_counts)
            self.stats[f'{name}_max_duplicate_count'] = duplicate_counts.max()
            
            # Analyze gene patterns in duplicates
            gene_col = 'Ref.Gene' if 'Ref.Gene' in df.columns else 'Gene.refGene'
            if gene_col in df.columns:
                genes_in_duplicates = duplicate_variants[gene_col].value_counts()
                self.stats[f'{name}_top_duplicated_genes'] = genes_in_duplicates.head(10).to_dict()
                
                # Log top genes with duplicates
                logging.info(f"{name}: Found {num_duplicates:,} duplicate variants affecting {len(duplicate_counts):,} unique positions")
                logging.info(f"{name}: Top 3 genes with duplicates: {dict(list(self.stats[f'{name}_top_duplicated_genes'].items())[:3])}")
        
        # Store basic statistics
        self.stats[f'{name}_total_variants'] = len(df)
        self.stats[f'{name}_unique_variants'] = len(variant_keys.unique())
        
        # Analyze chromosome distribution
        chr_dist = df['Chr'].value_counts().sort_index()
        self.stats[f'{name}_chr_distribution'] = chr_dist.to_dict()
        
        logging.info(f"{name} validation complete: {len(df):,} total variants, {len(variant_keys.unique()):,} unique variants")

    def read_files(self) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """
        Reads and validates input files.
        Handles column standardization and initial validation.
        """
        try:
            # Read InterVar file
            intervar_df = pd.read_csv(
                self.intervar_file, 
                sep='\t',
                low_memory=False
            )
            
            # Read MultiAnno file
            multianno_df = pd.read_csv(
                self.multianno_file,
                sep='\t',
                low_memory=False
            )
            
            # Standardize column names
            intervar_df.columns = [col.replace('#', '') for col in intervar_df.columns]
            multianno_df.columns = [col.replace('#', '') for col in multianno_df.columns]
            
            # Validate both dataframes
            self.validate_dataframe(intervar_df, 'InterVar')
            self.validate_dataframe(multianno_df, 'MultiAnno')
            
            logging.info(f"Successfully read input files. InterVar: {len(intervar_df):,} rows, MultiAnno: {len(multianno_df):,} rows")
            return intervar_df, multianno_df
            
        except Exception as e:
            logging.error(f"Error reading input files: {str(e)}")
            raise

    def validate_merge_results(self, merged_df: pd.DataFrame, intervar_df: pd.DataFrame) -> None:
        """
        Validates merge results and checks annotation completeness.
        """
        self.stats['total_merged_variants'] = len(merged_df)
        self.stats['variants_with_full_annotation'] = len(merged_df.dropna(subset=self.multianno_columns))
        self.stats['variants_with_partial_annotation'] = len(merged_df) - self.stats['variants_with_full_annotation']
        
        # Check for potential data loss
        if len(merged_df) != len(intervar_df):
            logging.warning(f"Merge resulted in different number of variants: InterVar={len(intervar_df):,}, Merged={len(merged_df):,}")
        
        # Check annotation completeness
        for col in self.multianno_columns:
            if col in merged_df.columns:
                missing = merged_df[col].isnull().sum()
                self.stats[f'missing_{col}'] = missing
                if missing > 0:
                    logging.info(f"Column {col}: {missing:,} variants lack annotation")

    def generate_summary_report(self) -> str:
        """
        Generates comprehensive summary report of the merge process.
        """
        report = [
            "=== Variant Merge Summary Report ===",
            "\nInput Statistics:",
            f"InterVar total variants: {self.stats.get('InterVar_total_variants', 0):,}",
            f"InterVar unique variants: {self.stats.get('InterVar_unique_variants', 0):,}",
            f"MultiAnno total variants: {self.stats.get('MultiAnno_total_variants', 0):,}",
            f"MultiAnno unique variants: {self.stats.get('MultiAnno_unique_variants', 0):,}",
            
            "\nDuplicate Analysis:",
            f"InterVar duplicates: {self.stats.get('InterVar_duplicates', 0):,}",
            f"Unique positions with duplicates: {self.stats.get('InterVar_unique_duplicated_variants', 0):,}",
            f"Maximum duplicates for a single variant: {self.stats.get('InterVar_max_duplicate_count', 0):,}",
            
            "\nChromosome Distribution (InterVar):"
        ]
        
        # Add chromosome distribution
        chr_dist = self.stats.get('InterVar_chr_distribution', {})
        for chr, count in sorted(chr_dist.items()):
            report.append(f"{chr}: {count:,} variants")
        
        # Add gene duplicate information if available
        if 'InterVar_top_duplicated_genes' in self.stats:
            report.extend([
                "\nTop Duplicated Genes:",
                "Gene: Number of variants"
            ])
            for gene, count in self.stats['InterVar_top_duplicated_genes'].items():
                report.append(f"{gene}: {count:,}")
        
        report.extend([
            "\nMerge Results:",
            f"Total merged variants: {self.stats.get('total_merged_variants', 0):,}",
            f"Variants with full annotation: {self.stats.get('variants_with_full_annotation', 0):,}",
            f"Variants with partial annotation: {self.stats.get('variants_with_partial_annotation', 0):,}"
        ])
        
        # Add annotation completeness details
        missing_annotations = [col for col in self.multianno_columns 
                             if f'missing_{col}' in self.stats and self.stats[f'missing_{col}'] > 0]
        if missing_annotations:
            report.append("\nMissing Annotations:")
            for col in missing_annotations:
                report.append(f"{col}: {self.stats[f'missing_{col}']:,} variants")
        
        return "\n".join(report)

    def merge_annotations(self) -> pd.DataFrame:
        """
        Performs the main annotation merge process.
        """
        try:
            # Read and validate input files
            intervar_df, multianno_df = self.read_files()
            
            # Create variant keys for matching
            intervar_df['variant_key'] = self.create_variant_key(intervar_df)
            multianno_df['variant_key'] = self.create_variant_key(multianno_df)
            
            # Select relevant columns from multianno
            available_columns = [col for col in self.multianno_columns if col in multianno_df.columns]
            multianno_subset = multianno_df[['variant_key'] + available_columns]
            
            # Merge dataframes
            merged_df = pd.merge(
                intervar_df,
                multianno_subset,
                on='variant_key',
                how='left'
            )
            
            # Drop the temporary variant key
            merged_df = merged_df.drop('variant_key', axis=1)
            
            # Fill missing values appropriately
            numeric_cols = ['REVEL_score', 'ClinPred_score', 'phyloP100way_vertebrate', 'GERP++_RS']
            for col in numeric_cols:
                if col in merged_df.columns:
                    merged_df[col] = merged_df[col].fillna(-1)
            
            categorical_cols = ['Polyphen2_HDIV_pred', 'MutationTaster_pred', 
                              'AlphaMissense_pred', 'Aloft_pred', 'PrimateAI_pred', 
                              'BayesDel_addAF_pred', 'CLNDN']
            for col in categorical_cols:
                if col in merged_df.columns:
                    merged_df[col] = merged_df[col].fillna('.')
            
            # Validate merge results
            self.validate_merge_results(merged_df, intervar_df)
            
            # Generate and log summary report
            summary_report = self.generate_summary_report()
            logging.info("\n" + summary_report)
            
            # Save summary report to file
            report_filename = "merge_summary_report.txt"
            with open(report_filename, 'w') as f:
                f.write(summary_report)
            logging.info(f"Summary report saved to {report_filename}")
            
            return merged_df
            
        except Exception as e:
            logging.error(f"Error during annotation merger: {str(e)}")
            raise

def main():
    """
    Main function to handle command line execution.
    """
    if len(sys.argv) != 4:
        print("Usage: python script.py intervarfile multiannofile output.tsv")
        sys.exit(1)
        
    intervar_file = sys.argv[1]
    multianno_file = sys.argv[2]
    output_file = sys.argv[3]
    
    try:
        # Create merger instance and process files
        merger = VariantMerger(intervar_file, multianno_file)
        merged_df = merger.merge_annotations()
        
        # Save output
        merged_df.to_csv(output_file, sep='\t', index=False)
        logging.info(f"Successfully wrote merged annotations to {output_file}")
        
    except Exception as e:
        logging.error(f"Error in main execution: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main()
