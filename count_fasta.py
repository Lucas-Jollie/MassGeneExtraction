#!/usr/bin/env python3

import argparse
import os
from Bio import SeqIO

def count_fasta_sequences(file_name):
    # Read the input FASTA file
    sequences = SeqIO.parse(file_name, "fasta")
    
    # Count the number of sequences
    count = sum(1 for _ in sequences)
    
    return count

def main():
    parser = argparse.ArgumentParser(description="Count the amount of FASTA sequences in a single file.")
    parser.add_argument('input_file', help='Input the FASTA file of which to count the number of sequences')
    
    args = parser.parse_args()
    # Count the number of sequences and print the result
    sequence_count = count_fasta_sequences(args.input_file)
    print(f"The file {args.input_file} contains {sequence_count} FASTA sequences.")
    
if __name__ == "__main__":
    main()