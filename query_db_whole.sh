#!/bin/bash

# Genes to loop over
genes=("vp1")
# loop through genes
for GENE in "${genes[@]}"; do
esearch -db nucleotide -query "($GENE[gene]) AND Enterovirus[organism]" | efetch -format gb > NCBI_whole_genomes_with_${GENE}.gb
done
