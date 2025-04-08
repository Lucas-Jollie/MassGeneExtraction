from Bio import SeqIO
from Bio.Seq import UndefinedSequenceError
import argparse


# Genes to loop over
genes = ["stx1A", 'stx1B', 'stx2A', 'stx2B']

def main():
    # Extracts all the genes in a given db
    parser = argparse.ArgumentParser(description="Extract annotated genes from a genbank database into a FASTA file.")
    parser.add_argument('input_file', help='Input Genbank database')
    parser.add_argument('output_file', help='Output FASTA file')
    
    args = parser.parse_args()

    input_file = args.input_file
    output_file = args.output_file

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