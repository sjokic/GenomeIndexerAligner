# -*- coding: utf-8 -*-
# @Time    : 11.10.20 20:35
# @Authors  : Chenxi-Nie / Guillaume-Thiry
# @E-mails    ：chenie@student.ethz.ch / gthiry@student.ethz.ch
# @FileName: bwa.py


## IMPORT

import numpy as np

## FUNCTIONS

def Burrows_Wheeler_Transform(sequence):
    '''
    Performs the BW transform using a quick sort
    :param sequence: genome sequence to index (string)
    :return: [bwt] the BW transform of the sequence
             [rotated_genome] the initial positions of the characters in the genome
    '''

    rotated_genome = [i for i in range(len(sequence))]

    # Sort the input sequencing based on quick sort
    def quick_sort(start, end, rotated_genome):

        # Trivial cases
        if end < start:
            return []

        if end == start:
            return [rotated_genome[end]]

        if end - start == 1:
            seperate_index = rotated_genome[start]
            sequence_1 = sequence[seperate_index:] + sequence[0:seperate_index]
            seperate_index = rotated_genome[end]
            sequence_2 = sequence[seperate_index:] + sequence[0:seperate_index]
            if sequence_1 < sequence_2:
                return rotated_genome[start: end + 1]
            else:
                return rotated_genome[start: end + 1][::-1]

        # General case

        seperate_index = rotated_genome[start] # Always use the first index to seperate rotated_genome into half
        seperate_rotation = sequence[seperate_index:] + sequence[0:seperate_index]
        seperate_index = start

        for i in range(start + 1, end + 1):
            ith_rotation = sequence[rotated_genome[i]:] + sequence[0:rotated_genome[i]]
            if ith_rotation < seperate_rotation:
                ith_index = rotated_genome[i]
                for j in range(i - 1, seperate_index - 1, -1):
                    rotated_genome[j + 1] = rotated_genome[j]
                rotated_genome[seperate_index] = ith_index
                seperate_index += 1

        rotated_genome_l = quick_sort(start, seperate_index - 1, rotated_genome)
        rotated_genome_r = quick_sort(seperate_index + 1, end,rotated_genome)

        return rotated_genome_l + [rotated_genome[seperate_index]] + rotated_genome_r


    rotated_genome = quick_sort(0, len(sequence) - 1, rotated_genome)
    #print(rotated_genome)

    bwt = ""
    for i in rotated_genome:
        if i == 0:
            bwt += sequence[-1]
        else:
            bwt += sequence[i - 1]
    return bwt, rotated_genome


def index_bwt(bwt):
    '''
    Given a BW transform, computes the Count and Occurence tables
    :param bwt: the BWT to consider (string)
    :return: [c] the count array of the different characters in the bwt
             [Occ] the occurence matrix of the bwt (one line per character)
    The order of the characters is : $ - A - C - G - T
    '''
    c = [0 for i in range(5)]
    for i in bwt:
        if i == "$":
            c[0] += 1
        elif i == "A":
            c[1] += 1
        elif i == "C":
            c[2] += 1
        elif i == "G":
            c[3] += 1
        elif i == "T":
            c[4] += 1

    num_A = c[1]
    c[0] = 0; c[1] = 1
    for i in [2, 3, 4]:
        temp = c[i]
        c[i] = c[i - 1] + num_A
        num_A = temp

    Occ = np.zeros((5, len(bwt)))
    for i in range(len(bwt)):
        if bwt[i] == "$":
            Occ[0][i : ] += 1
        elif bwt[i] == "A":
            Occ[1][i : ] += 1
        elif bwt[i] == "C":
            Occ[2][i : ] += 1
        elif bwt[i] == "G":
            Occ[3][i : ] += 1
        elif bwt[i] == "T":
            Occ[4][i : ] += 1

    return c, Occ


def search_in_bwt(bwt, sub_string):
    '''
    Uses the BWT of a sequence to find a given substring in the original sequence
    :param bwt: the BW transform (string)
    :param sub_string: the substring to find in the sequence
    :return: [sp, ep] left and right bound of the range in the rotated_genome where the substring is

    Reference : Algorithm and data structure of populational scale genomics lecture 02 P14
    '''

    c, Occ = index_bwt(bwt)
    p = len(sub_string)
    i = p - 1
    acgt_to_1234 = {"$":0, "A":1, "C":2, "G":3, "T":4}
    character = acgt_to_1234[sub_string[p - 1]]
    #print(c)
    #print(Occ)
    sp = c[character]
    ep = 0
    # quick fix not sure if works 100% correctly..

    if character == 4:
        ep = len(bwt) - 1
    else:
        ep = c[character + 1] - 1

    #ep = c[character + 1] - 1


    while(sp < ep and i >= 1):
        character = acgt_to_1234[sub_string[i - 1]]
        sp = int(c[character] + Occ[character][sp - 1])
        ep = int(c[character] + Occ[character][ep] - 1)
        i = i - 1
    return sp, ep


## EXAMPLES
'''
bwt, rotated_genome = Burrows_Wheeler_Transform("TAGAGAT$")
#print(bwt)
#print(rotated_genome)
sp, ep = search_in_bwt(bwt, "GAT")
print(sp)
print(ep)
for i in range(sp , ep + 1):
    print(rotated_genome[i])
'''
