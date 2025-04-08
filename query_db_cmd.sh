#!/bin/bash

# Check if the input file is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 input_file"
    exit 1
fi

input_file="$1"

# Check if the input file exists
if [ ! -f "$input_file" ]; then
    echo "File not found: $input_file"
    exit 1
fi

# Genes to loop over
# genes=("asta" "pet" "siga" "sat" "pic" "hlye")

while IFS= read -r GENE; do
    echo "Searching for gene: $GENE"
	esearch -db nucleotide -query "($GENE[gene]) AND Escherichia coli[organism] NOT (complete genome) NOT (contig) NOT (whole genome shotgun)" | efetch -format fasta > NCBI_${GENE}_sequences.fasta
done < "$input_file"

