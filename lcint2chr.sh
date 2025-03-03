#!/bin/bash

# Check if input file is provided
if [ $# -ne 1 ]; then
    echo "Usage: $0 <input_bed_file>"
    exit 1
fi

input_file=$1
output_file="${input_file%.bed}_chr.bed"

# Process the file: add 'chr' prefix to chromosome numbers
awk 'BEGIN{OFS="\t"} {
    # If first column is a number or X/Y/M, add "chr" prefix
    if ($1 ~ /^[0-9XYM]/) {
        $1 = "chr" $1
    }
    print
}' "$input_file" > "$output_file"

echo "Conversion complete! Output written to $output_file"
