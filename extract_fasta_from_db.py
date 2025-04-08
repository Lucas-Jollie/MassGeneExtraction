from Bio import SeqIO
from Bio.Seq import UndefinedSequenceError
import argparse
import os

# A function to go over each file in a directory to process them individually
# Currently I do not have a way to integrate all the files directly in case multiple
# files containing a gene are created. In the future they should be merged automatically

def process_directory(directory_path):
    
    # Loop through all files in the directory
    for root, dirs, files in os.walk(directory_path):
        for file in files:
            file_path = os.path.join(root, file)
            # file_path += str(n)
            create_fastas(file_path)
    

# Read gene names to loop over
# Now general because I mgiht use it for species as well
def read_input_files(file_of_interest):
    arg_list = open(file_of_interest,'r').read().splitlines()
    return arg_list

# Creates fasta files from the genes in the given gene list
def create_fastas(input_file):
    filename = input_file.split('.')
    filename = filename[0]
    
    # create the output file name
    output_file = filename+'_extracted'+'.fasta'
    # Open output file
    n = 0
    with open(input_file, "r") as input_handle, open(output_file, "w") as output_handle:
        # process each record in the database
        print(n)
        for record in SeqIO.parse(input_handle, "genbank"):
            # Go over all annotated features
            for feature in record.features:
                # Save features in gene list to FASTA
                if feature.type == "gene" and "gene" in feature.qualifiers:
                    gene_name = feature.qualifiers["gene"][0]
                    try:
                        gene_seq = feature.extract(record.seq)
                        output_handle.write(f">{gene_name}\n{gene_seq}\n")
                    except UndefinedSequenceError:
                        print(f"Sequence content is undefined for gene {gene_name} in {input_file}")


def main():
    # Ask and parse input
    parser = argparse.ArgumentParser(description="Extract annotated genes from a genbank database into\
    a FASTA file.")
    parser.add_argument('input_file', help='Input Genbank database file or directory containing the files')
    # parser.add_argument('output_file_handle', help='What to place at the start of the output file name for the created FASTA file')
    # parser.add_argument('gene_list_file', help='List of gene names to extract from the source')
    
    args = parser.parse_args()
    
    # Parse arguments to vars for readability
    # output_file_handle = args.output_file_handle
    # genes = read_input_files(args.gene_list_file)
    input_file = args.input_file
    
    # Check input for type and process accordingly
    if os.path.isfile(input_file):
        create_fastas(input_file)
    elif os.path.isdir(input_file):
        process_directory(input_file)
    else:
        print("Invalid input")

                        
if __name__ == "__main__":
    main()