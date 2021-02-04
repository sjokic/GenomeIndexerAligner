# -*- coding: utf-8 -*-
# @Time    : 23.10.20
# @Authors  : Chenxi N. / Guillaume T.
# @FileName: input.py

## IMPORT

import os
import time

## FUNCTIONS

def read_fa(file):
    '''
    Reads a FASTA file (.fa) containing a genome sequence
    :param file: full path of the .fa file
    :return: a string
    '''

    time_start = time.time()
    f = open(file)
    f.readline()
    genome = ""
    line = f.readline()
    while (line != ""):
        genome += line[:-1]
        line = f.readline()
    time_stop = time.time()
    print("Loading time : ", str(time_stop-time_start)[0:5], " s")
    return genome




def read_fq(file):
    '''
    Reads a FASTQ file (.fq) containing read sequences
    :param file: full path of the .fq file
    :return: a list of strings
    '''

    time_start = time.time()
    f = open(file)
    reads = []
    qnames = []
    quals = []
    b = True
    while(b):
        qname = f.readline()
        seq = f.readline()
        f.readline()
        qual = f.readline()
        if seq == "":
            b = False
        else:
            reads.append(seq[:-1])
            qnames.append(qname[1 : len(qname) - 3])
            quals.append(qual[:-1])
    time_stop = time.time()

    i = 0
    while qnames[0][i] != '-':
        i += 1

    print("Loading time : ", str(time_stop-time_start)[0:5], " s")
    return reads, qnames, qnames[0][0:i], quals

def read_sam(file, withPQline):
    '''
    Reads a Sam file (.sam)
    :param file: full path of the .fq file, boolean whether file contains PQ line
    :return: list of qname, flag, rname, pos, cigar
    '''

    time_start = time.time()
    f = open(file)
    qnames = []
    flags = []
    rname = []
    cigar = []
    pos = []
    print(f)
    f.readline()
    f.readline()
    if (withPQline):
        f.readline()

    while(True):
        line = f.readline()
        if line == "":
            break
        content = line.split('\t')
        #qnames.append(content[0])
        #flags.append(content[1])
        #rname.append(content[2])
        pos.append((content[3]))
        cigar.append(content[5])
    time_stop = time.time()

    print("Loading time : ", str(time_stop-time_start)[0:5], " s")
    return pos, cigar

## EXAMPLES

#x, y = read_sam("/home/adrian/PycharmProjects/BiomedicalEngineering/team09/data_small/output_tiny_30xCov.sam", True)

# Paths

#fa_path = os.getcwd() + "/data_small/genome.chr22.5K.fa"
#fq_path = os.getcwd() + "/data_small/output_tiny_30xCov1.fq"

# Loading

#genome = read_fa(fa_path)
#reads, qnames, sn = read_fq(fq_path)
#print(sn)
