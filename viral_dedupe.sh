#!/bin/bash

# set -x

# Check if the correct number of arguments is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <folders>"
    exit 1
fi

# Assign the arguments to variables
folder_files=$1 # File with bacterial species

# Initialize an array to hold the file contents
folders=()

# Read the file line by line and store each line in the array
while IFS= read -r gene; do
	# Checks if a new line was added in vain and ignores accordingly
	if [[ -n "$gene" ]]; then
		folders+=("$gene")
	fi
done < "$folder_files"

for folder in "${folders[@]}"; do
	echo $folder
	pwd
	cd "$folder/WholeGenomes"
	search_dir=$(pwd)
	echo "search dir $search_dir"
	for entry in "$search_dir"/*; do
		# echo "in search dirs $entry"
		outfile="${entry##$search_dir/}"
		outfile="deduped_$outfile"
		# exit 1
		bash dedupe.sh "in=$entry" "out=$outfile" "minidentity=95"
	done
	pwd
	cd "../.."
	pwd
done

# cd 'Norovirus'
# bash dedupe.sh

# cd 'Sapovirus'

# cd 'Enterovirus'

# cd 'HEV'

# cd 'Kobuvirus'

# cd 'Rotavirus'
