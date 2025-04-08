#!/bin/bash

# Hard coded file path untill all scripts are loaded into PATH
script_path='../../../../Python_Command_Tools/MassGeneExtraction/'
# script_path2='../../../../Python_Command_Tools/'

# Check if the correct number of arguments is provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <filename> <organism> <output_folder>"
    exit 1
fi

# Assign the arguments to variables
filename=$1 # File with bacterial species
organism=$2 # Organism to query for
output_folder=$3 # What directoty to save to (blank for current)

# Check and convert Windows line endings (CRLF) to Linux line endings (LF)
if file "$filename" | grep -q "CRLF"; then
    echo "Converting Windows line endings to Linux line endings..."
    sed -i 's/\r$//' "$filename"
fi

# Add a newline to the file to prevent file gene from being missed
echo "" >> "$filename"

# Create the output folder but silence the creation of it (-p)
mkdir -p $output_folder

# Initialize an array to hold the file contents
gene_list=()

# Read the file line by line and store each line in the array
while IFS= read -r gene; do
	# Checks if a new line was added in vain and ignores accordingly
	if [ -z "$gene" ]; then
		echo "Empty line found. Ignoring."
	# Adds actual txt to array
	else
		gene_list+=("$gene")
	fi
done < "$filename"

# Run the automated database search for each gene
for gene in "${gene_list[@]}"; do
	echo "Looking for whole genomes with $gene..."
	bash ${script_path}automated_query_db_whole_no_cmd.sh "$gene" "$organism" "$output_folder"
	echo "Extracting genes into FASTA format..."
	python3 ${script_path}extract_fasta_from_db.py "${output_folder}/NCBI_${gene}.gb"
	echo "Searching for individual gene uploads..."
	bash ${script_path}query_db_individual_genes.sh "$gene" "$organism" "$output_folder"
	echo "Extracting genes from individual querry into FASTA format..."
	python3 ${script_path}extract_fasta_from_db.py "${output_folder}/NCBI_${gene}_genes.gb"
	echo "Merging all genes into one file..."
	python3 ${script_path}combine_fasta_cmd_line.py "${output_folder}NCBI_${gene}_extracted.fasta" "${output_folder}NCBI_${gene}_genes_extracted.fasta" "${output_folder}NCBI_${gene}_merged.fasta"
done

# Create FASTA files from the databases you have just created
# python3 ${script_path}extract_fasta_from_db.py $output_folder

# Remove the (large) database files
echo "Removing database files..."
find "$output_folder" -type f -name "*.gb" -exec rm {} \;
echo "Removing bloated FASTA files..."
find "$output_folder" -type f -name "*extracted.fasta" -exec rm {} \;

# For each gene extract the sequences from the larger sets
for gene in "${gene_list[@]}"; do
	echo "Filtering FASTA files..."
	echo $gene
	filtered_file_name=$(grep -rl "$gene" "$output_folder" | head -n 1)
	echo $filtered_file_name
	python3 ${script_path}filter_fasta_cmd.py "$filtered_file_name" "${output_folder}NCBI_${gene}_filtered.fasta" "$gene" "True"
done

# remove the FASTA files containing all of the sequences to leave the extracted sequences
echo "Deleting merged files..."
find "$output_folder" -type f -name "*merged*" -exec rm {} \;


####TODO
# Use dedupe.sh to remove duplicates
# Not sure if it is wise to include dedupe.sh as it is not my code. It will make my use significantly easier.
