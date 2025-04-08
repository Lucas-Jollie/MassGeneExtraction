#!/bin/bash

# Hard coded file path untill all scripts are loaded into PATH
script_path='../../../../Python_Command_Tools/MassGeneExtraction/'
script_path2='../../../../Python_Command_Tools/'

# Check if the correct number of arguments is provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <filename> <organism> <output_folder>"
    exit 1
fi

# Assign the arguments to variables
filename=$1
organism=$2
output_folder=$3

mkdir $output_folder

# Initialize an array to hold the file contents
gene_list=()

# Read the file line by line and store each line in the array
while IFS= read -r gene; do
    gene_list+=("$gene")
done < "$filename"

echo $script_path
# Create FASTA files from the databases you have just created
python3 ${script_path}extract_fasta_from_db3.py $output_folder

exit 1

# Remove the (large) database files
find "$output_folder" -type f -name "*.gb" -exec rm {} \;

# For each gene extract the sequences from the larger sets
for gene in "${gene_list[@]}"; do
	filtered_file_name=$(grep -rl "$gene" "$output_folder" | head -n 1)
	echo $filtered_file_name
	python3 ${script_path2}filter_fasta_cmd.py "$filtered_file_name" "${output_folder}NCBI_${gene}_filtered.fasta" "$gene" "True"
done

# remove the FASTA files containing all of the sequences to leave the extracted sequences
find "$output_folder" -type f -name "*extracted*" -exec rm {} \;

## TODO
# Create something for individual querries
# Merge with output from individual querries
# Reformat individual querries
# Use dedupe.sh to remove duplicates
# 
