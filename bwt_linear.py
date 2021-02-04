# -*- coding: utf-8 -*-
# @Time    : 29.10.20 11:21
# @Author  : Guillaume T.
# @Updater : Chenxi N.
# @FileName: bwt_linear_time.py
# @E-mail    ：gthiry@student.ethz.ch / Chenie@student.ethz.ch

# Constructing Burrows Wheeler Transform
# in O(n) time

import numpy as np
import os
import time
from input import read_fa

from SA_linear import Suffix_Array

### CONSTRUCTION FUNCTIONS

def suffix_array_to_bwt(suffix_array, T):
    '''
    Transforms the Suffix Array of a sequence T into the corresponding BWT
    '''
    bwt = ""
    for i in range(len(T)):
        if suffix_array[i] > 0:
            bwt += T[suffix_array[i] - 1]
        else:
            bwt += "#"

    return bwt


convert_dict = {'#' : 1, 'A' : 2, 'C' : 3, 'G' : 4, 'N' : 5, 'T' : 6}

def bwt_linear(genome):
    '''
    Creates a BWT of the genome sequence in linear time
    '''

    # Initializing the input for the SA algorithm
    T = []
    n = len(genome)
    for i in range(n):
        T.append(convert_dict[genome[i]])
    T += [0,0,0]

    # Constructing SA in linear time
    SA = Suffix_Array(T,n,6)
    SA = [n] + SA

    # Converting SA to BWT in linear time
    bwt = suffix_array_to_bwt(SA, genome + '#')


    # Writing SA and BWT to files
    print("Writing bwt and suffix array to files")
    bwt_file = open("bwt.txt", "w")
    bwt_file.write(bwt)
    sa_file = open("suffix_array.txt", "w")
    sa_file.write(str(SA))

    return bwt, SA


### SEARCH FUNCTIONS


def index_bwt(bwt):
    c = [0 for i in range(6)]
    for i in bwt:
        if i == "$" or i == "#":
            c[0] += 1
        elif i == "A":
            c[1] += 1
        elif i == "C":
            c[2] += 1
        elif i == "G":
            c[3] += 1
        elif i == "N":
            c[4] += 1
        elif i == "T":
            c[5] += 1

    num_A = c[1]
    c[0] = 0; c[1] = 1
    for i in [2, 3, 4, 5]:
        temp = c[i]
        c[i] = c[i - 1] + num_A
        num_A = temp

    Occ = np.zeros((6, len(bwt)))
    start = time.time()
    for i in range(len(bwt)):
        if i > 0:
            Occ[0][i] = Occ[0][i - 1]
            Occ[1][i] = Occ[1][i - 1]
            Occ[2][i] = Occ[2][i - 1]
            Occ[3][i] = Occ[3][i - 1]
            Occ[4][i] = Occ[4][i - 1]
            Occ[5][i] = Occ[5][i - 1]

        if bwt[i] == "$" or bwt[i] == "#":
            Occ[0][i] += 1
        elif bwt[i] == "A":
            Occ[1][i] += 1
        elif bwt[i] == "C":
            Occ[2][i] += 1
        elif bwt[i] == "G":
            Occ[3][i] += 1
        elif bwt[i] == "N":
            Occ[4][i] += 1
        elif bwt[i] == "T":
            Occ[5][i] += 1
    return c, Occ


def search_in_bwt(bwt, sub_string, c, Occ):
    '''
        Uses the BWT of a sequence to find a given substring in the original sequence
        :param bwt: the BW transform (string)
        :param sub_string: the substring to find in the sequence
        :return: [sp, ep] left and right bound of the range in the rotated_genome where the substring is

        Reference : Algorithm and data structure of populational scale genomics lecture 02 P14
    '''

    p = len(sub_string)
    i = p - 1
    acgt_to_1234 = {"$":0, "A":1, "C":2, "G":3, "N":4,"T":5}
    character = acgt_to_1234[sub_string[p - 1]]

    sp = c[character]

    if character == 5:
        ep = len(bwt) - 1
    else:
        ep = c[character + 1] - 1


    while(sp <= ep and i >= 1):
        character = acgt_to_1234[sub_string[i - 1]]
        sp = int(c[character] + Occ[character][sp - 1])
        ep = int(c[character] + Occ[character][ep] - 1)
        i = i - 1
    return sp, ep


def search(bwt, sub_string, suffix_array, c, Occ):
    sp, ep = search_in_bwt(bwt, sub_string, c, Occ)

    hit = []
    for i in range(sp, ep + 1):
        hit.append(suffix_array[i])

    return hit


### EXAMPLES


# Step1 : load genome from fasta file
# genome = "TAGAGAT"

# # Step2 : build bwt and suffix array
# start = time.time()
# genome_bwt, genome_sa = bwt_linear(genome)
# end = time.time()
# print(str(end-start)[0:7], " s")

# Step3 : index the bwt before searching
# c, Occ = index_bwt(genome_bwt)

#sub_string = "GAT"
#hit = search(genome_bwt, sub_string, genome_sa, c, Occ)

# Do something with the hit
#print(hit)
