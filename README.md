# Genome indexer and aligner #

Genome indexer and aligner based on BWT for indexing (FM Index) and dynamic programming for alignment

Run main.py with the following arguments and in this particular order: `python main.py [name of directory containing .fa and .fq files] [name of genome .fa file] [name of .fq files]`

Example: `python main.py data_small genome.chr22.5K output_tiny_30xCov`

The output SAM file will be named "our_[name of .fq files].sam" and stored in the same directory as the .fa and .fq files.
