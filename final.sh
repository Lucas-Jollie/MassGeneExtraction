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
# echo "" >> "$filename"
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

len_of_genes=${#gene_list[@]}

if [ $len_of_genes == 1 ]; then
	endnote="FASTA file containing gene of interest is stored in ${output_folder}"
else
	endnote="FASTA files containing genes of interest are stored in ${output_folder}"
fi

echo $endnote
####TODO
# Use dedupe.sh to remove duplicates
# Not sure if it is wise to include dedupe.sh as it is not my code. It will make my use significantly easier.
