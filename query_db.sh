#!/bin/bash

# Genes to loop over
genes=("vp1")
# loop through genes
for GENE in "${genes[@]}"; do
esearch -db nucleotide -query "($GENE[gene]) AND Escherichia coli[organism] NOT (complete genome) NOT (contig) NOT (whole genome shotgun)" | efetch -format fasta > NCBI_${GENE}_sequences.fasta
done
