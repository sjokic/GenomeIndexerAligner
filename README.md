# Genome indexer and aligner #

Run main.py with the following arguments and in this particular order: `python main.py [name of directory containing .fa and .fq files] [name of genome .fa file] [name of .fq files]`

Example: `python main.py data_small genome.chr22.5K output_tiny_30xCov`

The output SAM file will be named "our_[name of .fq files].sam" and stored in the same directory as the .fa and .fq files.

### Important ###

For the chromosome 22 data located in the "data" directory, we have provided you with two text files named bwt_full.txt and sa_full.txt. These contain the Burrows–Wheeler transform and suffix array of chromosome 22 respectively. When running main.py you will be prompted to enter the names of these .txt files. Thus, when running main.py on choromsome 22, please enter "bwt_full" and "sa_full" respectively (without the .txt ending) if you wish to save time by skipping the computation of the BWT and suffix array. Otherwise enter anything else to continue and the BWT and suffix array will be computed instead.