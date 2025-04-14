#!/usr/bin/env python3
from Bio import SeqIO
import os
import argparse

def filter_fasta(input_file, output_file, keyword, flag):
    # Read the input FASTA file
    sequences = SeqIO.parse(input_file, "fasta")
    
    if flag:
        # Filter sequences that do contain the keyword in their title
        filtered_sequences = (seq for seq in sequences if keyword.lower() in seq.description.lower())
        # filtered_sequences = (seq for seq in sequences if keyword.lower() in seq.description)
    elif flag==False:
        # Filter sequences that do not contain the keyword in their title
        print("False flag triggered")
        filtered_sequences = (seq for seq in sequences if keyword.lower() not in seq.description.lower())
    else:
        print('Invalid input')
    
    # Write the filtered sequences to the output FASTA file
    SeqIO.write(filtered_sequences, output_file, "fasta")

# print(os.getcwd())

# Prompt the user for input and output file names
# input_file = input("Enter the input FASTA file name: ")
# output_file = input("Enter the output FASTA file name: ")
# keyword = input("Enter the keyword to filter from fasta names: ")

# input_file = "Bacteria/Food/Campylobacter/ena_coding_Ccoli_flaB_with_hyoilei.fasta" 
# output_file = "Bacteria/Food/Campylobacter/ena_coding_Ccoli_flaB_filtered.fasta"
# keyword = "hyoilei"



def main():
    parser = argparse.ArgumentParser(description="Filter sequences from FASTA file using a keyword found in the headers.")
    parser.add_argument('input_file', help='Input FASTA file')
    parser.add_argument('output_file', help='Output FASTA file')
    parser.add_argument('keyword', help='Headers containing this word are left out.')
    parser.add_argument('keep_key', nargs='?', type=bool, default=False, help='Keep the keyword?')

    args = parser.parse_args()

    filter_fasta(args.input_file, args.output_file, args.keyword, args.keep_key)
    print(f"Filtered FASTA file created: {args.output_file}")

if __name__ == "__main__":
    main()