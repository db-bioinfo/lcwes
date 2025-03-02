#!/usr/bin/env python3
import pandas as pd
import json
from datetime import datetime
from variant_prioritizer import VariantPrioritizer
import sys
import os
import traceback
import numpy as np

def process_variant_file(input_file: str, output_file: str) -> None:
    """Process variant file and generate prioritized output"""
    try:
        # Read the input file
        print(f"Reading input file: {input_file}")
        df = pd.read_csv(input_file, sep='\t', low_memory=False)
        print(f"Successfully loaded {len(df)} variants")
        
        # Initialize prioritizer
        prioritizer = VariantPrioritizer()
        
        # Process variants
        print("Processing variants...")
        df_prioritized = prioritizer.process_variants(df)
        
        # Take only the top 25000 variants
        df_top_variants = df_prioritized.head(25000)
        print(f"Limiting output to top 25,000 variants")
        
        # Generate stats filename
        stats_file = output_file.rsplit('.', 1)[0] + '_stats.json'
        
        # Save prioritized variants
        print(f"Saving prioritized variants to: {output_file}")
        df_top_variants.to_csv(output_file, sep='\t', index=False)
        
        # Generate and save statistics
        print("Generating statistics...")
        stats = generate_statistics(df_top_variants)
        print(f"Saving statistics to: {stats_file}")
        with open(stats_file, 'w') as f:
            json.dump(stats, f, indent=2)
        
        # Print summary
        print("\nProcessing complete!")
        print("\nVariant Classification Summary:")
        for tier, count in stats['classification_summary'].items():
            print(f"{tier}: {count['count']} variants ({count['percentage']:.2f}%)")
        
        print(f"\nTotal variants processed: {stats['total_variants']}")
        print(f"Mean priority score: {stats['score_statistics']['mean']:.3f}")
        print(f"Median priority score: {stats['score_statistics']['median']:.3f}")
        print(f"\nOriginal variants count: {len(df_prioritized)}")
        print(f"Top variants saved: {len(df_top_variants)}")
            
    except Exception as e:
        print(f"Error processing file: {str(e)}")
        print("\nDetailed error traceback:")
        print(traceback.format_exc())
        sys.exit(1)

def generate_statistics(df: pd.DataFrame) -> dict:
    """Generate summary statistics for the processed variants"""
    try:
        total_variants = len(df)
        
        # Classification summary
        classification_counts = df['classification'].value_counts()
        classification_summary = {}
        for tier in classification_counts.index:
            count = classification_counts[tier]
            classification_summary[tier] = {
                'count': int(count),
                'percentage': (count/total_variants) * 100
            }
        
        # Score statistics
        stats = {
            'total_variants': total_variants,
            'classification_summary': classification_summary,
            'score_statistics': {
                'mean': float(df['final_score'].mean()),
                'median': float(df['final_score'].median()),
                'min': float(df['final_score'].min()),
                'max': float(df['final_score'].max()),
                'quartiles': {
                    'q1': float(df['final_score'].quantile(0.25)),
                    'q2': float(df['final_score'].quantile(0.50)),
                    'q3': float(df['final_score'].quantile(0.75))
                }
            },
            'high_impact_variants': {
                'total': int(len(df[df['final_score'] >= 0.8])),
                'top_genes': df[df['final_score'] >= 0.8]['Ref.Gene'].value_counts().head(10).to_dict()
            },
            'timestamp': datetime.now().isoformat()
        }
        
        return stats
        
    except Exception as e:
        print(f"Error generating statistics: {str(e)}")
        return {
            'total_variants': len(df),
            'error': str(e),
            'timestamp': datetime.now().isoformat()
        }

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python prioritize_variants.py <input_file> <output_file>")
        print("Example: python prioritize_variants.py input.txt output.tsv")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    # Check if input file exists
    if not os.path.exists(input_file):
        print(f"Error: Input file '{input_file}' does not exist.")
        sys.exit(1)
    
    # Check if output directory exists
    output_dir = os.path.dirname(output_file)
    if output_dir and not os.path.exists(output_dir):
        print(f"Creating output directory: {output_dir}")
        try:
            os.makedirs(output_dir)
        except Exception as e:
            print(f"Error creating output directory: {str(e)}")
            sys.exit(1)
    
    # Process variants
    process_variant_file(input_file, output_file)
