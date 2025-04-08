README Version 0.1.1
Script versions 0.2
Author: Lucas Jollie, MSc
University of Amsterdam, MBMFS

There are two options to run:
1. Run for individual genes
2. Download the entire genomes of organisms and extract sequences of interest from there

Running both and combining them may give you more hits in total. Future versions will be updated to use similar
filenames to create ease of use.

query_db_whole.sh
This script will download entire genbank databases from which target genes can be extracted using
extract_fasta_from_db.py
The script unfortunately for now has to be run by hardcoding the organism and genes of interest into the query.
Alternatively run the query through the command line using environmental variables. However, the idea is this can be run for a set of organisms and genes to automatically download massive amounts of data easily.

extract_fasta_from_db.py
This Python script will extract the gene(s) of interest from the db you created using query_db_whole.sh
Future updates will allow either use of the command line or an input file with genes of interest. For now,
manually update the gene names by hardcoding them into the script.
