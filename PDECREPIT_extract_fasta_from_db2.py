from Bio import SeqIO
from Bio.Seq import UndefinedSequenceError
import argparse


# Read gene names to loop over
def read_input_files(file_of_interest):
    arg_list = open(file_of_interest,'r').read().splitlines()
    return arg_list


def main():
    parser = argparse.ArgumentParser(description="Extract annotated genes from a genbank database into\
    a FASTA file.")
    parser.add_argument('input_file', help='Input Genbank database')
    parser.add_argument('output_file_handle', help='What to place at the start of the output file name \
    for the created FASTA file')
    parser.add_argument('gene_list_file', help='List of gene names to extract from the source')
    
    args = parser.parse_args()

    input_file = args.input_file
    output_file_handle = args.output_file_handle
    genes = read_input_files(args.gene_list_file)

    for gene in genes:
        output_file = output_file_handle+'_'+gene+'.fasta'
        with open(input_file, "r") as input_handle, open(output_file, "w") as output_handle:
            for record in SeqIO.parse(input_handle, "genbank"):
                for feature in record.features:
                    if feature.type == "gene" and "gene" in feature.qualifiers:
                        gene_name = feature.qualifiers["gene"][0]
                        try:
                            gene_seq = feature.extract(record.seq)
                            output_handle.write(f">{gene_name}\n{gene_seq}\n")
                        except UndefinedSequenceError:
                            print(f"Sequence content is undefined for gene {gene_name} in {input_file}")

                        
if __name__ == "__main__":
    main()