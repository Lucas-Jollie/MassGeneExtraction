#!/bin/bash

# Check if the input file is provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 gene organism output_folder"
    exit 1
fi

gene="$1"
organism="$2"
output_folder="$3"

echo "Starting individual gene query..."
echo "Searching for gene: $gene"
esearch -db nucleotide -query "(${gene}[gene]) AND ${organism}[organism] NOT (complete genome) NOT (contig) NOT (whole genome shotgun)" | efetch -format gb > ${output_folder}NCBI_${gene}_genes.gb

