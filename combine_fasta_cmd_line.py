#!/usr/bin/env python3

import os
import argparse

def combine_fasta_files(input_files, output_file):
    with open(output_file, 'w') as outfile:
        for file in input_files:
            with open(file, 'r') as infile:
                for line in infile:
                    outfile.write(line)

def main():
    parser = argparse.ArgumentParser(description="Combine multiple FASTA files into one.")
    parser.add_argument('input_files', nargs='+', help='List of input FASTA files')
    parser.add_argument('output_file', help='Output FASTA file')

    args = parser.parse_args()

    combine_fasta_files(args.input_files, args.output_file)
    print(f"Combined FASTA file created: {args.output_file}")

if __name__ == "__main__":
    main()