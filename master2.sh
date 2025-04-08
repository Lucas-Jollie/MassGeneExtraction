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
filename=$1 # File with bacterial species
organism=$2 # Organism to query for
output_folder=$3 # What directoty to save to (blank for current)

# Create the output folder but silence the creation of it (-p)
mkdir -p $output_folder

# Initialize an array to hold the file contents
gene_list=()

# Read the file line by line and store each line in the array
while IFS= read -r gene; do
    gene_list+=("$gene")
done < "$filename"

# Run the automated database search for each gene
for gene in "${gene_list[@]}"; do
	echo "Looking for whole genomes with $gene..."
	# bash ${script_path}automated_query_db_whole_no_cmd.sh "$gene" "$organism" "$output_folder"
	echo "Extracting genes into FASTA format..."
	python3 ${script_path}extract_fasta_from_db.py "${output_folder}/NCBI_${gene}.gb"
	echo "Searching for individual gene uploads..."
	bash ${script_path}query_db_individual_genes.sh "$gene" "$organism" "$output_folder"
	echo "Merging all genes into one file..."
	python3 ${script_path2}combine_fasta_cmd_line.py "${output_folder}NCBI_${gene}_extracted.fasta" "${output_folder}NCBI_${gene}_sequences.fasta" "${output_folder}NCBI_${gene}_merged.fasta"
done

# Create FASTA files from the databases you have just created
# python3 ${script_path}extract_fasta_from_db.py $output_folder

# Remove the (large) database files
echo "Removing database files..."
find "$output_folder" -type f -name "*.gb" -exec rm {} \;
echo "Removing bloated FASTA files..."
find "$output_folder" -type f -name "*extracted.fasta" -exec rm {} \;
echo "Removing individual sequences..."
find "$output_folder" -type f -name "*sequences.fasta" -exec rm {} \;
# do I want to add the merger here so I can filter at the same time?

# For each gene extract the sequences from the larger sets
for gene in "${gene_list[@]}"; do
	echo "Filtering FASTA files..."
	echo $gene
	filtered_file_name=$(grep -rl "$gene" "$output_folder" | head -n 1)
	echo $filtered_file_name
	python3 ${script_path2}filter_fasta_cmd.py "$filtered_file_name" "${output_folder}NCBI_${gene}_filtered.fasta" "$gene" "True"
done

# remove the FASTA files containing all of the sequences to leave the extracted sequences
# find "$output_folder" -type f -name "*extracted*" -exec rm {} \;



# Reformat individual querries
# Use dedupe.sh to remove duplicates
# Not sure if it is wise to include dedupe.sh as it is not my code. It will make my use significantly easier.
# I don't know how to tackle issues where I have an entire operon yet.
# E.g. stx1a is often combined with stx1b so I need to get that out. Happens for (some) other genes too.