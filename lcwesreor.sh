#!/bin/bash

# Extract columns for clinical reporting: Chr Start End Ref Alt Transcript:HGVS VariantType
# Usage: ./clinical_hgvs_with_type.sh input.tsv > output.tsv

awk -F'\t' '{
    # Extract the transcript ID and c. notation
    split($3, parts, ",")
    first_anno = parts[1]
    if (match(first_anno, /NM_[0-9]+/)) {
        nm_id = substr(first_anno, RSTART, RLENGTH)
        if (match(first_anno, /c\.[^:]+/)) {
            c_notation = substr(first_anno, RSTART, RLENGTH)
            transcript_hgvs = nm_id ":" c_notation
        } else {
            transcript_hgvs = first_anno
        }
    } else {
        transcript_hgvs = first_anno
    }
    
    # Include the variant type (column 2) in the output
    print $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" transcript_hgvs "\t" $2
}' "$1"
