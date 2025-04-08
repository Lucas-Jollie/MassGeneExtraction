# For each gene extract the sequences from the larger sets
for gene in "${gene_list[@]}"; do
	filtered_file_name=$(grep -rl "$gene" "$output_folder" | head -n 1)
	echo $filtered_file_name
	python3 ${script_path2}filter_fasta_cmd.py "$filtered_file_name" "${output_folder}NCBI_${gene}_filtered.fasta" "$gene" "True"
done

# remove the FASTA files containing all of the sequences to leave the extracted sequences
find "$output_folder" -type f -name "*extracted*" -exec rm {} \;

# Merge with output from individual querries
for gene in "${gene_list[@]}"; do
	newfile="(${output_folder}NCBI_${gene}_merged.fasta)"
	grep -l "$gene_name" "$directory"/*.fasta | xargs cat > "$newfile"

done

# Reformat individual querries
# Use dedupe.sh to remove duplicates
# Not sure if it is wise to include dedupe.sh as it is not my code. It will make my use significantly easier.
# I don't know how to tackle issues where I have an entire operon yet.
# E.g. stx1a is often combined with stx1b so I need to get that out. Happens for (some) other genes too.