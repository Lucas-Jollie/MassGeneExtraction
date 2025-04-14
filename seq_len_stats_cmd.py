from Bio import SeqIO
import argparse

def calculate_average_length(fasta_file):
    total_length = 0
    sequence_count = 0
    min_len = 0
    max_len = 0
    
    for record in SeqIO.parse(fasta_file, "fasta"):
        rec_len = len(record.seq)
        total_length += rec_len
        sequence_count += 1
        if min_len == 0 or min_len > rec_len:
            min_len = rec_len
        if max_len < rec_len:
            max_len = rec_len
    
    average_length = total_length / sequence_count if sequence_count > 0 else 0
    return average_length, min_len, max_len, sequence_count

def main():
    parser = argparse.ArgumentParser(description="Determine min, max and avg sequence lengths.")
    parser.add_argument('input_file', help='Input FASTA file')

    args = parser.parse_args()

    # combine_fasta_files(args.input_files, args.output_file)
    
    average_length, min_length, max_length, num_seqs = calculate_average_length(args.input_file)
    print(f"Total amount of sequences: {num_seqs} sequences")
    print(f"Minimun sequence length: {min_length} bp")
    print(f"Maximum sequence length: {max_length} bp")
    print(f"Average sequence length: {average_length} bp")

    # print(f"Combined FASTA file created: {args.output_file}")

if __name__ == "__main__":
    main()

# fasta_file = input("Enter the FASTA file of which you want the average sequence length: ")
# average_length, min_length, max_length = calculate_average_length(fasta_file)
# print(d"Minimun sequence length: {min_length} bp")
# print(d"Maximum sequence length: {max_length} bp")
# print(f"Average sequence length: {average_length} bp")
