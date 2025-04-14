#!/bin/bash

# set -x

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

# Remove any trailing newlines
sed -i -e :a -e '/^\n*$/{$d;N;};/\n$/ba' "$filename"

# Add a single newline at the end
echo "" >> "$filename"

# Create the output folder but silence the creation of it (-p)
mkdir -p $output_folder

# Check and convert Windows line endings (CRLF) to Linux line endings (LF)
if file "$filename" | grep -q "CRLF"; then
    echo "Converting Windows line endings to Linux line endings..."
    sed -i 's/\r$//' "$filename"
fi

if [[ "${output_folder: -1}" != "/" ]]; then
	output_folder="${output_folder}/"
fi

# Initialize an array to hold the file contents
gene_list=()

# Read the file line by line and store each line in the array
while IFS= read -r gene; do
	# Checks if a new line was added in vain and ignores accordingly
	if [[ -n "$gene" ]]; then
		gene_list+=("$gene")
	fi
done < "$filename"

# Run the automated database search for each gene
for gene in "${gene_list[@]}"; do
	echo "Looking for whole genomes with $gene..."
	bash automated_query_db_whole_no_cmd.sh "$gene" "$organism" "$output_folder"
	echo "Extracting genes from whole genome databases into FASTA format..."
	extract_fasta_from_db.py "${output_folder}/NCBI_${gene}.gb"
	echo "Searching for individual gene uploads..."
	bash query_db_individual_genes.sh "$gene" "$organism" "$output_folder"
	echo "Extracting genes from individual query into FASTA format..."
	extract_fasta_from_db.py "${output_folder}/NCBI_${gene}_genes.gb"
	echo "Merging all genes into one file..."
	combine_fasta_cmd_line.py "${output_folder}NCBI_${gene}_extracted.fasta" "${output_folder}NCBI_${gene}_genes_extracted.fasta" "${output_folder}NCBI_${gene}_merged.fasta"
	echo "Files merged succesfully!"
	echo "End of loop for $gene"
done

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
	filter_fasta_cmd.py "$filtered_file_name" "${output_folder}NCBI_${gene}_filtered.fasta" "$gene" "True"
done

# remove the FASTA files containing all of the sequences to leave the extracted sequences
echo "Deleting merged files..."
find "$output_folder" -type f -name "*merged*" -exec rm {} \;

len_of_genes=${#gene_list[@]}

if [ $len_of_genes == 1 ]; then
	endnote="FASTA file containing gene of interest for $organism is stored in ${output_folder}"
else
	endnote="FASTA files containing genes of interest for $organism are stored in ${output_folder}"
fi

echo $endnote
####TODO
# Use dedupe.sh to remove duplicates
# Not sure if it is wise to include dedupe.sh as it is not my code. It will make my use significantly easier.
