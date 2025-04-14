#!/usr/bin/env python3
# Import BioPython and stats functions
from Bio import SeqIO
import argparse
from seq_len_stats_cmd import calculate_average_length
# using the 

def filter_by_length(fasta, results, criterium, multiplier_min=4, multiplier_max=3):
    multiplier_max = float(multiplier_max)
    multiplier_min = float(multiplier_min)
    
    with open(results, 'w') as output_handle:
        for record in SeqIO.parse(fasta, "fasta"):
            if criterium/multiplier_min < len(record.seq) < multiplier_max*criterium:
                SeqIO.write(record, output_handle, "fasta")

def main():

    parser = argparse.ArgumentParser(description="Filter sequences based on average length.")
    parser.add_argument('input_file', help='Input FASTA file')
    parser.add_argument('output_file', help='Output FASTA file')
    parser.add_argument('min_bound_multiplier', nargs='?', const=4, help='Input FASTA file')
    parser.add_argument('max_bound_multiplier', nargs='?', const=3, help='Input FASTA file')

    args = parser.parse_args()
    input_seqs = args.input_file
    average_length = calculate_average_length(input_seqs)[0]
    filter_by_length(input_seqs, args.output_file, average_length,
    args.min_bound_multiplier, args.max_bound_multiplier)


if __name__ == "__main__":
    main()
