#!/usr/bin/env python3
"""
High-Performance Variant Annotation Merger

This script merges annotations from a clnsnp file into a mutect file,
matching variants by chromosome, position, and alleles with enhanced
capabilities for handling complex variants and different representations.

Usage:
  python lcwesmergeClnSnp2Prior.py <mutect_file> <clnsnp_file> <output_file> [--threads N]

Author: Clinical Bioinformatics Team
"""

import sys
import csv
import re
import logging
import argparse
import os
import time
from collections import defaultdict
import functools
import concurrent.futures
import multiprocessing
import itertools
import array
import heapq
from multiprocessing import shared_memory
import pickle
import io
import threading
from contextlib import contextmanager

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[logging.StreamHandler()]
)
logger = logging.getLogger(__name__)

# Lock for thread-safe progress updates
progress_lock = threading.Lock()

# Global counters for progress tracking
TOTAL_VARIANTS = 0
PROCESSED_VARIANTS = 0
MATCHED_VARIANTS = 0
START_TIME = 0

# Global shared data dictionary (used by threads)
SHARED_DATA = {}

# Function to memoize expensive operations
def memoize(func):
    """Cache results of function calls to avoid repeated computation."""
    cache = {}
    @functools.wraps(func)
    def wrapper(*args):
        if args in cache:
            return cache[args]
        result = func(*args)
        cache[args] = result
        return result
    return wrapper

@memoize
def left_normalize_indel(ref, alt, pos):
    """
    Normalize indel representation by left-trimming common prefixes
    and right-trimming common suffixes.
    
    Args:
        ref (str): Reference allele
        alt (str): Alternate allele
        pos (int): Position
        
    Returns:
        tuple: (normalized_ref, normalized_alt, normalized_pos)
    """
    # Handle special cases
    if ref == alt:
        return ref, alt, pos
    
    if ref == '-':
        return '', alt, pos
    
    if alt == '-':
        return ref, '', pos
    
    # Left-trim common prefix
    i = 0
    while i < min(len(ref), len(alt)) and ref[i] == alt[i]:
        i += 1
    
    if i > 0:
        ref = ref[i:]
        alt = alt[i:]
        pos += i
    
    # Right-trim common suffix
    j = 0
    while j < min(len(ref), len(alt)) and ref[-(j+1)] == alt[-(j+1)]:
        j += 1
    
    if j > 0:
        ref = ref[:-j]
        alt = alt[:-j]
    
    # If we've trimmed everything, represent as empty strings
    if not ref:
        ref = ''
    if not alt:
        alt = ''
    
    return ref, alt, pos

@memoize
def parse_rs_id(rs_id):
    """Extract numeric RS ID from string with 'rs' prefix."""
    if not rs_id or rs_id == '.':
        return None
    match = re.search(r'rs(\d+)', rs_id)
    return match.group(1) if match else None

def get_gene_aliases(gene):
    """Get all possible aliases for a gene."""
    # Simply return the gene itself as aliases functionality has been removed
    return [gene]

def create_lookup_indexes(clnsnp_variants):
    """
    Create optimized lookup indexes for variant matching.
    
    Args:
        clnsnp_variants (list): List of variants from clnsnp file
        
    Returns:
        dict: Dictionary of indexes for different matching strategies
    """
    # Pre-allocate dictionaries with estimated capacity
    estimated_variants = len(clnsnp_variants)
    
    logger.info(f"Creating optimized indexes for {estimated_variants:,} variants...")
    
    # Primary indexes for fast lookups
    exact_matches = {}
    normalized_matches = {}
    rs_id_matches = {}
    pos_gene_matches = {}  # For position±1 matching
    indel_matches = {}     # For indel position±3 matching
    gene_proximity = {}    # For gene in proximity matching
    hgvs_matches = {}      # For HGVS notation matching
    multi_allelic = {}     # For expanded multi-allelic variants
    
    # Process each variant
    for variant_idx, variant in enumerate(clnsnp_variants):
        # Extract key fields
        chrom = variant['CHROM']
        try:
            pos = int(variant['POS'])
        except ValueError:
            pos = 0
        ref = variant['REF']
        alt = variant['ALT']
        gene = variant.get('ANN[0].GENE', '')
        rs_id = variant.get('RS', '').strip()
        
        # Index 1: Exact matching by chrom:pos:ref:alt
        exact_key = f"{chrom}:{pos}:{ref}:{alt}"
        exact_matches[exact_key] = variant_idx
        
        # Check for multi-allelic variants
        is_multi_allelic = ',' in alt
        if is_multi_allelic:
            # Expand multi-allelic variants
            alt_alleles = alt.split(',')
            
            # Check if VAF is also comma-separated
            vaf_values = []
            if 'VAF' in variant and ',' in variant['VAF']:
                vaf_values = variant['VAF'].split(',')
            
            # Create normalized representation for each allele
            for i, alt_allele in enumerate(alt_alleles):
                # Make a shallow copy of relevant data for the expanded variant
                new_variant = variant.copy()
                new_variant['ALT'] = alt_allele
                
                # Update VAF if available
                if vaf_values and i < len(vaf_values):
                    new_variant['VAF'] = vaf_values[i]
                
                # Get normalized representation
                if ref == '-':
                    norm_ref = ''
                else:
                    norm_ref = ref
                    
                if alt_allele == '-':
                    norm_alt = ''
                else:
                    norm_alt = alt_allele
                
                norm_ref, norm_alt, norm_pos = left_normalize_indel(norm_ref, norm_alt, pos)
                norm_key = f"{chrom}:{norm_pos}:{norm_ref}:{norm_alt}"
                
                # Add to multi-allelic index
                multi_allelic[norm_key] = new_variant
        
        # Proceed with other indexes if not multi-allelic
        if not is_multi_allelic:
            # Index 2: RS ID matching
            if rs_id and rs_id != '.':
                rs_id_matches[rs_id] = variant_idx
            
            # Index 3: Normalized representation
            if ref == '-':
                norm_ref = ''
            else:
                norm_ref = ref
                
            if alt == '-':
                norm_alt = ''
            else:
                norm_alt = alt
                
            norm_ref, norm_alt, norm_pos = left_normalize_indel(norm_ref, norm_alt, pos)
            norm_key = f"{chrom}:{norm_pos}:{norm_ref}:{norm_alt}"
            normalized_matches[norm_key] = variant_idx
        
        # Gene-based indexes (apply to all variants)
        if gene:
            # Index 4: Position and gene matching for tolerance ±1
            for offset in range(-1, 2):  # -1, 0, 1
                pos_gene_key = f"{chrom}:{pos+offset}:{gene}"
                # Don't overwrite if the key already exists (preserve smaller offset)
                if pos_gene_key not in pos_gene_matches:
                    pos_gene_matches[pos_gene_key] = (abs(offset), variant_idx)
                elif abs(offset) < pos_gene_matches[pos_gene_key][0]:
                    pos_gene_matches[pos_gene_key] = (abs(offset), variant_idx)
            
            # Index 5: Indel position matching for larger tolerance (±3)
            if len(ref) > 1 or len(alt) > 1 or ref == '-' or alt == '-':
                for offset in range(-3, 4):  # -3 to +3
                    if offset in [-1, 0, 1]:  # Already covered above
                        continue
                    indel_key = f"{chrom}:{pos+offset}:{gene}:indel"
                    # Don't overwrite if the key already exists (preserve smaller offset)
                    if indel_key not in indel_matches:
                        indel_matches[indel_key] = (abs(offset), variant_idx)
                    elif abs(offset) < indel_matches[indel_key][0]:
                        indel_matches[indel_key] = (abs(offset), variant_idx)
            
            # Index 6: Gene proximity (±5 bp)
            for offset in range(-5, 6):
                if offset in range(-3, 4):  # Already covered above
                    continue
                prox_key = f"{chrom}:{pos+offset}:{gene}"
                # Don't overwrite if the key already exists (preserve smaller offset)
                if prox_key not in gene_proximity:
                    gene_proximity[prox_key] = (abs(offset), variant_idx)
                elif abs(offset) < gene_proximity[prox_key][0]:
                    gene_proximity[prox_key] = (abs(offset), variant_idx)
            
            # Index 7: HGVS notation
            hgvs_c = variant.get('ANN[0].HGVS_C', '')
            if hgvs_c:
                hgvs_key = f"{gene}:{hgvs_c}"
                hgvs_matches[hgvs_key] = variant_idx
    
    # Finalize the position-based indexes
    # Convert tuple (offset, variant_idx) to just variant_idx in the final maps
    pos_gene_matches_final = {k: v[1] for k, v in pos_gene_matches.items()}
    indel_matches_final = {k: v[1] for k, v in indel_matches.items()}
    gene_proximity_final = {k: v[1] for k, v in gene_proximity.items()}
    
    # Create the final indexes dictionary
    indexes = {
        'exact': exact_matches,
        'normalized': normalized_matches,
        'rs_id': rs_id_matches,
        'pos_gene': pos_gene_matches_final,
        'indel': indel_matches_final,
        'proximity': gene_proximity_final,
        'hgvs': hgvs_matches,
        'multi_allelic': multi_allelic
    }
    
    logger.info(f"Created optimized indexes: "
               f"exact={len(exact_matches):,}, "
               f"normalized={len(normalized_matches):,}, "
               f"rs_id={len(rs_id_matches):,}, "
               f"multi_allelic={len(multi_allelic):,}")
    
    return indexes, clnsnp_variants

def find_best_match(mutect_variant, indexes, clnsnp_variants):
    """
    Find the best matching variant in clnsnp for the given mutect variant using
    optimized lookup indexes.
    
    Args:
        mutect_variant (dict): Variant from mutect file
        indexes (dict): Dictionary of lookup indexes
        clnsnp_variants (list): Original clnsnp variants
    
    Returns:
        dict: Best matching variant from clnsnp, or None if no match found
    """
    # Extract key fields for matching
    chrom = mutect_variant['Chr']
    try:
        pos = int(mutect_variant['Start'])
    except ValueError:
        pos = 0
    ref = mutect_variant['Ref']
    alt = mutect_variant['Alt']
    gene = mutect_variant['Ref.Gene']
    rs_id = parse_rs_id(mutect_variant.get('avsnp151', None))
    
    # Normalize ref and alt for matching
    if ref == '-':
        ref = ''
    if alt == '-':
        alt = ''
    
    # Get normalized representation for comparison
    norm_ref, norm_alt, norm_pos = left_normalize_indel(ref, alt, pos)
    
    # Strategy 1: Direct position and allele match (fastest)
    exact_key = f"{chrom}:{pos}:{ref}:{alt}"
    if exact_key in indexes['exact']:
        variant_idx = indexes['exact'][exact_key]
        return clnsnp_variants[variant_idx]
    
    # Strategy 2: Match normalized variants
    norm_key = f"{chrom}:{norm_pos}:{norm_ref}:{norm_alt}"
    if norm_key in indexes['normalized']:
        variant_idx = indexes['normalized'][norm_key]
        return clnsnp_variants[variant_idx]
    
    # Strategy 3: Check multi-allelic variants
    if norm_key in indexes['multi_allelic']:
        return indexes['multi_allelic'][norm_key]
    
    # Strategy 4: Match by RS ID
    if rs_id and rs_id in indexes['rs_id']:
        variant_idx = indexes['rs_id'][rs_id]
        return clnsnp_variants[variant_idx]
    
    # Strategy 5: Match by position with tolerance of ±1
    for offset in [-1, 0, 1]:
        pos_gene_key = f"{chrom}:{pos+offset}:{gene}"
        if pos_gene_key in indexes['pos_gene']:
            variant_idx = indexes['pos_gene'][pos_gene_key]
            return clnsnp_variants[variant_idx]
    
    # Strategy 6: Match indels with larger position tolerance (±3)
    if len(ref) > 1 or len(alt) > 1 or ref == '-' or alt == '-':
        for offset in range(-3, 4):
            if offset in [-1, 0, 1]:  # Already checked
                continue
            indel_key = f"{chrom}:{pos+offset}:{gene}:indel"
            if indel_key in indexes['indel']:
                variant_idx = indexes['indel'][indel_key]
                return clnsnp_variants[variant_idx]
    
    # Strategy 7: Match by gene in proximity (±5 bp)
    for offset in range(-5, 6):
        if offset in range(-3, 4):  # Already checked
            continue
        prox_key = f"{chrom}:{pos+offset}:{gene}"
        if prox_key in indexes['proximity']:
            variant_idx = indexes['proximity'][prox_key]
            return clnsnp_variants[variant_idx]
    
    # Strategy 8: Match by HGVS notation
    otherinfo = mutect_variant.get('Otherinfo', '')
    hgvs_match = re.search(r'c\.([A-Za-z0-9_>.+-]+)', otherinfo)
    if hgvs_match and gene:
        hgvs_c = f"c.{hgvs_match.group(1)}"
        hgvs_key = f"{gene}:{hgvs_c}"
        if hgvs_key in indexes['hgvs']:
            variant_idx = indexes['hgvs'][hgvs_key]
            return clnsnp_variants[variant_idx]
    
    # No match found
    return None

def process_chunk(chunk_data):
    """
    Process a chunk of mutect variants in parallel.
    
    Args:
        chunk_data (tuple): (chunk_id, variants, columns_to_add)
        
    Returns:
        tuple: (processed_variants, matched_count)
    """
    global SHARED_DATA, PROCESSED_VARIANTS, MATCHED_VARIANTS, TOTAL_VARIANTS, START_TIME
    
    chunk_id, variants, columns_to_add = chunk_data
    indexes = SHARED_DATA['indexes']
    clnsnp_variants = SHARED_DATA['clnsnp_variants']
    
    results = []
    matched = 0
    
    # Process each variant in the chunk
    for variant in variants:
        # Find the best match
        best_match = find_best_match(variant, indexes, clnsnp_variants)
        
        # Add columns from the match or use '.' for missing data
        for column in columns_to_add:
            if best_match and column in best_match:
                variant[column] = best_match[column].strip('\r')
                if column == columns_to_add[0]:  # Count match based on first column
                    matched += 1
            else:
                variant[column] = "."
        
        results.append(variant)
    
    # Update global progress counters
    with progress_lock:
        PROCESSED_VARIANTS += len(variants)
        MATCHED_VARIANTS += matched
        
        # Report progress periodically but don't slow down processing
        current_time = time.time()
        elapsed = current_time - START_TIME
        
        # Calculate stats for progress reporting
        progress_pct = PROCESSED_VARIANTS / TOTAL_VARIANTS * 100 if TOTAL_VARIANTS > 0 else 0
        match_rate = MATCHED_VARIANTS / PROCESSED_VARIANTS * 100 if PROCESSED_VARIANTS > 0 else 0
        variants_per_sec = PROCESSED_VARIANTS / elapsed if elapsed > 0 else 0
        remaining_variants = TOTAL_VARIANTS - PROCESSED_VARIANTS
        remaining_time = remaining_variants / variants_per_sec if variants_per_sec > 0 else 0
        
        # Report only from one thread occasionally to avoid log spam
        if chunk_id % 5 == 0:
            logger.info(
                f"Progress: {PROCESSED_VARIANTS:,}/{TOTAL_VARIANTS:,} variants "
                f"({progress_pct:.1f}%), "
                f"Match rate: {match_rate:.1f}%, "
                f"Speed: {variants_per_sec:.1f} variants/sec, "
                f"Est. remaining: {int(remaining_time/60)}m {int(remaining_time%60)}s"
            )
    
    return results, matched

def read_file_with_progress(filename, description, process_func):
    """Read a file with progress reporting."""
    result = None
    line_count = 0
    start_time = time.time()
    last_report_time = start_time
    report_interval = 2  # seconds
    
    try:
        # First, try to quickly count lines for progress reporting
        logger.info(f"Counting lines in {filename}...")
        with open(filename, 'r') as f:
            for _ in f:
                line_count += 1
        logger.info(f"File contains {line_count:,} lines")
        
        # Now process the file
        logger.info(f"Reading {description} from {filename}...")
        with open(filename, 'r') as f:
            processed_lines = 0
            result = process_func(f)
            
            # Report completion
            elapsed = time.time() - start_time
            rate = line_count / elapsed if elapsed > 0 else 0
            logger.info(f"Completed reading {description}: {line_count:,} lines in {elapsed:.1f} seconds ({rate:.1f} lines/sec)")
    
    except Exception as e:
        logger.error(f"Error reading {filename}: {str(e)}")
        raise
    
    return result

def read_clnsnp_file(filename):
    """Read and parse clnsnp file."""
    def process(file_obj):
        reader = csv.DictReader(file_obj, delimiter='\t')
        return [{k.strip('\r'): v.strip('\r') for k, v in row.items()} for row in reader]
    
    return read_file_with_progress(filename, "clnsnp variants", process)

def read_mutect_file(filename):
    """Read and parse mutect file."""
    def process(file_obj):
        reader = csv.DictReader(file_obj, delimiter='\t')
        return list(reader), reader.fieldnames
    
    return read_file_with_progress(filename, "mutect variants", process)

def main(mutect_file, clnsnp_file, output_file, num_threads=None, chunk_size=None, use_threads=True):
    """
    Main function to process and merge variant annotations.
    
    Args:
        mutect_file (str): Path to mutect file
        clnsnp_file (str): Path to clnsnp file
        output_file (str): Path to output file
        num_threads (int, optional): Number of threads to use
        chunk_size (int, optional): Size of chunks for processing
        use_threads (bool): Whether to use threads instead of processes
    """
    global SHARED_DATA, TOTAL_VARIANTS, PROCESSED_VARIANTS, MATCHED_VARIANTS, START_TIME
    
    # Initialize global counters
    TOTAL_VARIANTS = 0
    PROCESSED_VARIANTS = 0
    MATCHED_VARIANTS = 0
    
    # Determine optimal thread count
    if num_threads is None:
        cpu_count = multiprocessing.cpu_count()
        # For large files, use more threads than CPUs since it's IO-bound
        if cpu_count >= 16:
            num_threads = cpu_count * 2
        else:
            num_threads = max(4, cpu_count)
    
    logger.info(f"Starting variant annotation merger using {num_threads} threads")
    
    # Columns to add from clnsnp
    columns_to_add = [
        "ANN[0].FEATUREID", 
        "ANN[0].EFFECT", 
        "ANN[0].HGVS_C", 
        "ANN[0].HGVS_P", 
        "Filter", 
        "CLNHGVS", 
        "VAF", 
        "RS"
    ]
    
    # Step 1: Read clnsnp file
    clnsnp_variants = read_clnsnp_file(clnsnp_file)
    logger.info(f"Loaded {len(clnsnp_variants):,} variants from clnsnp file")
    
    # Step 2: Create optimized lookup indexes
    indexes, clnsnp_variants = create_lookup_indexes(clnsnp_variants)
    
    # Step 3: Read mutect file
    mutect_variants, fieldnames = read_mutect_file(mutect_file)
    logger.info(f"Loaded {len(mutect_variants):,} variants from mutect file")
    
    # Set total variants for progress tracking
    TOTAL_VARIANTS = len(mutect_variants)
    
    # Step 4: Determine optimal chunk size
    if chunk_size is None:
        # For very large files, use larger chunks
        if len(mutect_variants) > 100000:
            # Balance between parallelism and overhead
            chunk_size = max(1000, len(mutect_variants) // (num_threads * 5))
        else:
            chunk_size = max(100, len(mutect_variants) // (num_threads * 10))
    
    # Split variants into chunks for processing
    chunks = []
    for i in range(0, len(mutect_variants), chunk_size):
        chunk = mutect_variants[i:i+chunk_size]
        chunks.append((i//chunk_size, chunk, columns_to_add))
    
    logger.info(f"Processing in {len(chunks)} chunks (chunk size: {chunk_size}) using {num_threads} threads")
    
    # Step 5: Share data with worker threads/processes
    SHARED_DATA['indexes'] = indexes
    SHARED_DATA['clnsnp_variants'] = clnsnp_variants
    
    # Step 6: Process chunks in parallel
    START_TIME = time.time()
    all_results = []
    total_matched = 0
    
    # Choose between ThreadPoolExecutor (faster, lower overhead) and ProcessPoolExecutor
    executor_class = concurrent.futures.ThreadPoolExecutor if use_threads else concurrent.futures.ProcessPoolExecutor
    
    with executor_class(max_workers=num_threads) as executor:
        futures = [executor.submit(process_chunk, chunk_data) for chunk_data in chunks]
        
        # Collect results as they complete
        for future in concurrent.futures.as_completed(futures):
            try:
                results, matched = future.result()
                all_results.extend(results)
                total_matched += matched
            except Exception as e:
                logger.error(f"Error processing chunk: {str(e)}")
                raise
    
    # Calculate final statistics
    processing_time = time.time() - START_TIME
    variants_per_second = PROCESSED_VARIANTS / processing_time if processing_time > 0 else 0
    
    logger.info(f"Parallel processing completed in {processing_time:.2f} seconds ({variants_per_second:.1f} variants/second)")
    logger.info(f"Matching complete: {MATCHED_VARIANTS:,} matched, {TOTAL_VARIANTS-MATCHED_VARIANTS:,} unmatched")
    logger.info(f"Match rate: {MATCHED_VARIANTS/TOTAL_VARIANTS*100:.1f}%")
    
    # Step 7: Write results to output file
    logger.info(f"Writing output to: {output_file}")
    write_start = time.time()
    
    with open(output_file, 'w', newline='') as outfile:
        # Create new header with added columns
        all_fieldnames = fieldnames + columns_to_add
        
        writer = csv.DictWriter(outfile, fieldnames=all_fieldnames, delimiter='\t')
        writer.writeheader()
        
        # Use larger batches for writing to improve performance
        write_batch_size = 10000
        total_written = 0
        
        for i in range(0, len(all_results), write_batch_size):
            batch = all_results[i:i+write_batch_size]
            for row in batch:
                writer.writerow(row)
            
            total_written += len(batch)
            if total_written % 100000 == 0:
                logger.info(f"Writing output: {total_written:,}/{len(all_results):,} variants written ({total_written/len(all_results)*100:.1f}%)")
    
    write_time = time.time() - write_start
    logger.info(f"Writing completed in {write_time:.2f} seconds")
    logger.info(f"Merged annotations written to {output_file}")

def run_from_command_line():
    """Parse command line arguments and run the script."""
    parser = argparse.ArgumentParser(description='High-Performance Variant Annotation Merger')
    parser.add_argument('mutect_file', help='Input mutect file')
    parser.add_argument('clnsnp_file', help='Input clnsnp file')
    parser.add_argument('output_file', help='Output merged file')
    parser.add_argument('--threads', type=int, help='Number of threads to use (default: auto-tuned based on system)')
    parser.add_argument('--chunk-size', type=int, help='Size of chunks for parallel processing (default: auto-tuned)')
    parser.add_argument('--processes', action='store_true', help='Use processes instead of threads (slower but more robust)')
    parser.add_argument('--verbose', action='store_true', help='Enable verbose logging')
    
    args = parser.parse_args()
    
    # Set logging level
    if args.verbose:
        logger.setLevel(logging.DEBUG)
    
    # Start timer for full execution
    start_time = time.time()
    
    # Run the main function
    main(
        args.mutect_file,
        args.clnsnp_file,
        args.output_file,
        args.threads,
        args.chunk_size,
        not args.processes
    )
    
    # Report total execution time
    total_time = time.time() - start_time
    hours = int(total_time // 3600)
    minutes = int((total_time % 3600) // 60)
    seconds = total_time % 60
    
    if hours > 0:
        time_str = f"{hours}h {minutes}m {seconds:.1f}s"
    elif minutes > 0:
        time_str = f"{minutes}m {seconds:.1f}s"
    else:
        time_str = f"{seconds:.2f} seconds"
    
    logger.info(f"Total execution time: {time_str}")

if __name__ == "__main__":
    run_from_command_line()
