# Overview
I downloaded 5,000 RSV fasta sequences+their metadata from GISAID. Sample collection dates span from December 2011 to April 2023.

In this project I analyze both the metadata and the sequence data

## Project structure
- **blast_results/**: Contains the blast alignment results for our RSV sequences
- **gisaid_data/**: Contains starting fasta files and metadata.
- **reference_fasta_files/**: Contains two fasta files used as references for my blast alignments.
- **blast_class.py**: file containing classes used to discern blast alignments results.
- **graph_rsv.py**: file containing various functions used to graph my results in jupyter notebook.
- **metadata_rsv.py**: file containing functions used to extract and discern information from the downloaded gisaid metadata.
- **rsv_overview.ipynb**: jupyter notebook which runs through all of the above.
