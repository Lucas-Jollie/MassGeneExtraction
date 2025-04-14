#!/bin/bash

# Hard coded file path untill all scripts are loaded into PATH
script_path='../../../../Python_Command_Tools/MassGeneExtraction/'

# Check if the correct number of arguments is provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <filename> <organism> <output_folder>"
    exit 1
fi

# Assign the arguments to variables
filename=$1 # File with genes
organism=$2 # Organism to query for
output_folder=$3 # What directoty to save to (blank for current)

# Check and convert Windows line endings (CRLF) to Linux line endings (LF)
# if file "$organism_files" | grep -q "CRLF"; then
    # echo "Converting Windows line endings to Linux line endings..."
    # sed -i 's/\r$//' "$organism_files"
# fi

# Add a newline to the file to prevent file gene from being missed
# echo "" >> "$organism_files"

echo "Moving to DAEC..."
cd 'DAEC'
bash master2.sh "$filename" "$organism" "$output_folder"

echo "Moving to EAEC..."
cd '../EAEC/'
bash master2.sh "$filename" "$organism" "$output_folder"

echo "Moving to EIEC..."
cd '../EIEC/'
bash master2.sh "$filename" "$organism" "$output_folder"

echo "Moving to EPEC..."
cd '../EPEC/'
bash master2.sh "$filename" "$organism" "$output_folder"

echo "Moving to ETEC..."
cd '../ETEC/'
bash master2.sh "$filename" "$organism" "$output_folder"

echo "Completed all E. coli variants"
