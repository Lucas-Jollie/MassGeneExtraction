#!/bin/bash

# Hardcoded paths for now
script_path='../../../../Python_Command_Tools/MassGeneExtraction/'

# Check if the correct number of arguments is provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <filename> <organism> <output_folder>"
    exit 1
fi

# Assign the arguments to variables
filename=$1 # File with bacterial species
organism=$2 # Organism to query for
output_folder=$3 # What directoty to save to (blank for current)

bash ${script_path}master.sh "$filename" "$organism" "$output_folder"

# Initialize an array to hold the file contents
gene_list=()

# Read the file line by line and store each line in the array
while IFS= read -r gene; do
    gene_list+=("$gene")
done < "$filename"

for gene in "${gene_list[@]}"; do
	queryfile="${output_folder}NCBI_${gene}_filtered.fasta"
	echo "deduping FASTA files..."
	bash dedupe.sh "$queryfile" "${output_folder}NCBI_${gene}_deduped.fasta" "minidentity=95"