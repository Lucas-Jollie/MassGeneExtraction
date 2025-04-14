**README Version 0.2.0
Script versions 1.0
Author:** Lucas Jollie, MSc
University of Amsterdam, MBMFS

# NCBI gene and genomic region downloader and extractor
This is a simple command line tool I developed to download and extract genes or genomic regions of interest
into FASTA files from the NCBI database. The script first looks for whole genomes annotated with the genes
of interest and stores them in a database file before extracting all annotated genes.
Next, other uploads of the genes of interest or CDS annotated with the genes of interest are downloaded into
a database format and extracted into FASTA format.
FASTA files are then merged and filtered to contain only the genes or genomic regions of interest.
Do **not** interfere with temporary files that are created during the running of the script or you will brick
the process. Database files and individual files are removed at the end anyway.

# System requirements:
The script is made to run on UNIX systems, but has only been tested on Ubuntu 22.04 as of yet.
This script requires Python3 and BioPython to be installed.
IMPORTANT: Make sure you have added the directory of scripts to your PATH, else they will not run.

# Installation
Pull the repository into your preferred location of choice. Make sure to add the folder to your PATH variable
if you have not done so already.

# Usage
* *$ master.sh [genes.txt] [organism name] [Output folder]* *
Currently the script only accepts text files with genes separated by newlines. The script is designed to work
with Unix-style line endings, but also accepts Windows-style line endings. 


## Individual scripts
### automated_query_db_whole_no_cmd.sh
This script will download entire genbank databases from which target genes can be extracted using
extract_fasta_from_db.py
The script unfortunately for now has to be run by hardcoding the organism and genes of interest into the query.
Alternatively run the query through the command line using environmental variables. However, the idea is this can be run for a set of organisms and genes to automatically download massive amounts of data easily.
Usage:
* *$ automated_query_db_whole_no_cmd.sh [genes.txt] [target organism] [output folder]* *

### query_db_individual_genes.sh
Runs exactly the same as the first script, but ignores whole genomes. Data is extracted in the same way as
for whole genomes.
Usage:
* *$ query_db_individual_genes.sh [genes.txt] [target organism] [output folder]* *

### extract_fasta_from_db.py
This Python script will extract the gene(s) of interest from the db you created using query_db_whole.sh
Future updates will allow either use of the command line or an input file with genes of interest. For now,
manually update the gene names by hardcoding them into the script.
Usage:
* *$ extract_fasta_from_db.sh [genbank_database.gb]* *

### combine_fasta_cmd_line.py
A python script that takes command line input to merge two or more FASTA files into one.
Usage:
* *$ python combine_fasta_cmd_line.py [input file(s)] [output file]* *

### filter_fasta_cmd.py
A python script that takes command line input of a FASTA file. The script can either save or delete the input.
Default in * *master.sh* * is to keep the input.

# Future perspectives and other options
## Other uses
### Navigating directories
This script can be combined or reiterated by the user with another script to run for multiple organisms
and directories for example. If each directory contains the same gene input file (e.g. genes.txt) a simple
script that navigates directories can read the input file and download the genes per species.

### Combining with other tools such as BBtools
I have combined this script and scripts from BBtools such as * *dedupe.sh* * to remove duplicates and cluster
output into one file, create antisense strands and even design DNA hybridisation probes. All credit to the
original authors of those tools.


## Future perspectives
- I might add the directory navigating option myself.
- I might also add functionality to delete duplicate entries, as you will now have to do so yourself manually or using other tools.
- Perhaps delete the files per gene query to stop the disk space from becoming to full of redundant files when downloading many genes
- 