# Authors: Stefan Jokic, Adrian Taubner

import os
import sys
import gzip
import time

from align import getCigarFromIndex, reverse_complement
from bwt_linear import bwt_linear, index_bwt, search
from input import read_fa, read_fq
import numpy as np

from writeOutput import writeToFile
from getFlags import getFlags


def run(path, readFile, genomeFile):
    read1 = path + readFile + "1.fq"
    read2 = path + readFile + "2.fq"
    # Get files
    if(not os.path.exists(path + genomeFile + ".fa")):
        gz = gzip.open(path + genomeFile + ".fa" + ".gz", 'rb')
        fq = open(path + genomeFile + ".fa", 'wb')
        fq.write(gz.read())

    if(not os.path.exists(read1)):
        gz = gzip.open(read1 + ".gz", 'rb')
        fq = open(read1, 'wb')
        fq.write(gz.read())

    if(not os.path.exists(read2)):
        gz = gzip.open(read2 + ".gz", 'rb')
        fq = open(read2, 'wb')
        fq.write(gz.read())

    print("Loading {} ...".format(genomeFile + ".fa"))
    genome = read_fa(path + genomeFile + ".fa")
    print("Loading {} ...".format(readFile + "1.fq"))
    reads1, qnamesRead1, sequenceName, quals1  = read_fq(read1)
    print("Loading {} ...".format(readFile + "2.fq"))
    reads2, qnamesRead2, sequenceName, quals2  = read_fq(read2)

    # initialize variables
    Read1_isFwd = True
    Read2_isFwd = True

    print("")
    bwt_txt = input("Is BWT of given genome already stored in a .txt file? If so specify its file name without the .txt ending.\n Otherwise enter anything to continue. \n")
    sa_txt = input("Is suffix array of given genome already stored in a .txt file? If so specify its file name without the .txt ending.\n Otherwise enter anything to continue. \n")
    print("")


    if((os.path.exists(bwt_txt + ".txt")) and (os.path.exists(sa_txt + ".txt"))):
        print("Loading BWT and suffix array from .txt files ...")
        time_start = time.time()
        bwt_file = open(bwt_txt + ".txt", 'r')
        genome_bwt = bwt_file.read()
        sa_file = open(sa_txt + ".txt", 'r')
        genome_sa = list(map(int, list(sa_file.read()[1:-1].split(", "))))
        time_stop = time.time()
        print("Finished loading BWT and suffix array from .txt files in " + str(time_stop-time_start)[0:5] + " s.\n")
    else:
        print("Constructing BWT and suffix array of genome ...")
        time_start = time.time()
        genome_bwt, genome_sa = bwt_linear(genome)
        time_stop = time.time()
        print("Finished constructing BWT and suffix array of genome in " + str(time_stop-time_start)[0:5] + " s.\n")


    c, Occ = index_bwt(genome_bwt)
    min_len_seed = 5
    outputPath = path + "our_" + readFile + ".sam"
    if(os.path.exists(outputPath)):
        inp = ''
        while(inp != 'y' or inp != 'n'):
            inp = input("WARNING: Output .sam file already exists! Are you sure you want to continue (File will be overwritten) ? (y/n)\n")
            if(inp == 'y'):
                os.remove(outputPath)
                break
            elif (inp == 'n'):
                sys.exit("Aborted.")
            else:
                print("Please enter 'y' to continue or 'n' to abort.")

    file = open(outputPath, 'a')
    header = True

    for i in range(0,len(reads1)):
        if(i == 1):
            header = False
        print("Processing read pair #" + str(i+1) + ", out of " + str(len(reads1)))
        read1_fwd = reads1[i]
        read1_rvcmpl = reverse_complement(read1_fwd)

        read2_fwd = reads2[i]
        read2_rvcmpl = reverse_complement(read2_fwd)

        seedLength1 = len(read1_fwd)
        seedLength2 = len(read1_fwd)

        hits1_fwd = search(genome_bwt, read1_fwd[0:seedLength1], genome_sa, c, Occ)
        hits1_rvcmpl = search(genome_bwt, read1_rvcmpl[0:seedLength1], genome_sa, c, Occ)

        hits2_fwd = search(genome_bwt, read2_fwd[0:seedLength2], genome_sa, c, Occ)
        hits2_rvcmpl = search(genome_bwt, read2_rvcmpl[0:seedLength2], genome_sa, c, Occ)


        while seedLength1 > min_len_seed and len(hits1_fwd) == 0 and len(hits1_rvcmpl) == 0:
            seedLength1 = int(seedLength1 / 2)

            if (seedLength1 < 2*min_len_seed):
                seedLength1 = min_len_seed
            for j in range(int(len(read1_fwd) / seedLength1)):
                hits1_fwd = search(genome_bwt, read1_fwd[seedLength1 * j:seedLength1 * (j + 1)], genome_sa, c, Occ)
                hits1_rvcmpl = search(genome_bwt, read1_rvcmpl[seedLength1 * j:seedLength1 * (j + 1)], genome_sa, c, Occ)
                if len(hits1_fwd) != 0 or len(hits1_rvcmpl) != 0:
                    break

        while seedLength2 > min_len_seed and len(hits2_fwd) == 0 and len(hits2_rvcmpl) == 0:
            seedLength2 = int(seedLength2 / 2)

            if (seedLength2 < 2*min_len_seed):
                seedLength2 = min_len_seed
            for j in range(int(len(read2_fwd) / seedLength2)):
                hits2_fwd = search(genome_bwt, read2_fwd[seedLength2 * j:seedLength2 * (j + 1)], genome_sa, c, Occ)
                hits2_rvcmpl = search(genome_bwt, read2_rvcmpl[seedLength2 * j:seedLength2 * (j + 1)], genome_sa, c, Occ)
                if len(hits2_fwd) != 0 or len(hits2_rvcmpl) != 0:
                    break


        minCostRead1 = sys.maxsize
        minCigarRead1 = ''
        minPosRead1 = 0

        minCostRead2 = sys.maxsize
        minCigarRead2 = ''
        minPosRead2 = 0

        isMapped1 = False
        isMapped2 = False

        if(len(hits1_fwd) or len(hits1_rvcmpl)):
            isMapped1 = True

        if(len(hits2_fwd) or len(hits2_rvcmpl)):
            isMapped2 = True

        # forward alignment
        for hit in hits1_fwd:
            leftIndex, rightIndex = 0, len(genome)
            leftIndex = max(hit - len(read1_fwd), leftIndex)
            rightIndex = min(hit + len(read2_rvcmpl) + 1, rightIndex)
            cigar, cost, index = getCigarFromIndex(read1_fwd, genome[leftIndex : rightIndex])
            if cost < minCostRead1:
                minCigarRead1 = cigar
                minPosRead1 = leftIndex + index + 1 # because MAPQ index starts at 1
                minCostRead1 = cost
        for hit in hits2_fwd:
            leftIndex, rightIndex = 0, len(genome)
            leftIndex = max(hit - len(read2_fwd), leftIndex)
            rightIndex = min(hit + len(read2_rvcmpl) + 1, rightIndex)
            cigar, cost, index = getCigarFromIndex(read2_fwd, genome[leftIndex : rightIndex])
            if cost < minCostRead2:
                minCigarRead2 = cigar
                minPosRead2 = leftIndex + index + 1 # because MAPQ index starts at 1
                minCostRead2 = cost

        # reverse complement alignment
        for hit in hits1_rvcmpl:
            leftIndex, rightIndex = 0, len(genome)
            leftIndex = max(hit - len(read1_rvcmpl), leftIndex)
            rightIndex = min(hit + len(read2_rvcmpl) + 1, rightIndex)
            cigar, cost, index = getCigarFromIndex(read1_rvcmpl, genome[leftIndex : rightIndex])
            if cost < minCostRead1:
                minCigarRead1 = cigar
                minPosRead1 = leftIndex + index + 1 # because MAPQ index starts at 1
                minCostRead1 = cost
                Read1_isFwd = False
        for hit in hits2_rvcmpl:
            leftIndex, rightIndex = 0, len(genome)
            leftIndex = max(hit - len(read2_rvcmpl), leftIndex)
            rightIndex = min(hit + len(read2_rvcmpl) + 1, rightIndex)
            cigar, cost, index = getCigarFromIndex(read2_rvcmpl, genome[leftIndex : rightIndex])
            if cost < minCostRead2:
                minCigarRead2 = cigar
                minPosRead2 = leftIndex + index + 1 # because MAPQ index starts at 1
                minCostRead2 = cost
                Read2_isFwd = False

        flag1, flag2 = getFlags(Read1_isFwd, Read2_isFwd, isMapped1, isMapped2)

        alignmentsRead1 = [qnamesRead1[i], flag1, minPosRead1, minCigarRead1, reads1[i], Read1_isFwd, quals1[i]]
        alignmentsRead2 = [qnamesRead2[i], flag2, minPosRead2, minCigarRead2, reads2[i], Read2_isFwd, quals2[i]]
        outputPath = path + "our_" + readFile + ".sam"
        writeToFile(file, header, sequenceName, len(genome), alignmentsRead1, alignmentsRead2)

    file.close()

# Paths
if((len(sys.argv) < 4) or (len(sys.argv) > 4)):
    sys.exit("ERROR: Incorrect number of arguments.\n Please specify the following arugments and in this particular order: Name of directory, name of genome fa file and name of fq files.\n Example: python main.py data_small genome.chr22.5K output_tiny_30xCov")

if((not os.path.exists(sys.argv[1])) or ((not os.path.exists(sys.argv[1] + "/" + sys.argv[2] + ".fa")) and (not os.path.exists(sys.argv[1] + "/" + sys.argv[2] + ".fa.gz"))) or (((not os.path.exists(sys.argv[1] + "/" + sys.argv[3] + "1.fq")) or (not os.path.exists(sys.argv[1] + "/" + sys.argv[3] + "2.fq"))) and ((not os.path.exists(sys.argv[1] + "/" + sys.argv[3] + "1.fq.gz")) or (not os.path.exists(sys.argv[1] + "/" + sys.argv[3] + "2.fq.gz"))))):
    sys.exit("ERROR: Could not find specified files.\n Please specify the following arugments and in this particular order: Name of directory, name of genome fa file and name of fq files.\n Example: python main.py data_small genome.chr22.5K output_tiny_30xCov")

fq_path = sys.argv[3]
fa_path = sys.argv[2]
folder = os.getcwd() + "/" + sys.argv[1] + "/"

'''
fq_path = "output_tiny_30xCov"
fa_path = "genome.chr22.5K"
folder = os.getcwd() + "/" + "data_small" + "/"
'''

run(folder, fq_path, fa_path)
