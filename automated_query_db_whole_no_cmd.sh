#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <gene> <organism> <output_folder>"
    exit 1
fi

# Assign the arguments to variables
gene=$1
organism=$2
output_folder=$3

#Initialize an array to hold the file contents
# gene_list=()

#Read the file line by line and store each line in the array
# while IFS= read -r gene; do
    # gene_list+=("$gene")
# done < "$filename"

# Execute the esearch and efetch commands for a given filename
echo "Processing gene: $gene"
echo "Organism: $organism"
query="(${gene}[gene] AND ${organism}[organism])"
# echo $query
esearch -db nucleotide -query "$query" | efetch -format gb > ${output_folder}NCBI_${gene}.gb

# Execute the esearch command
# esearch -db nucleotide -query "(${gene}[gene]) AND ${organism}[organism]" | efetch -format gb > "${output_folder}NCBI_whole_genomes_with_${gene}.gb"

