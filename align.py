# Authors: Stefan Jokic, Adrian Taubner

import numpy as np

mismatch = 2
gap = 2.5
match = -1

def reverse_complement(sequence):
    '''
    Generates the reverse complement of a sequence.
    :param sequence: input sequence.
    :return: reverse complement of input sequence.
    '''
    comp = {'A' : 'T', 'C' : 'G', 'G' : 'C', 'T' : 'A'}
    rev_comp = ''

    for i in range(0, len(sequence)):
        rev_comp = comp[sequence[i]] + rev_comp

    return rev_comp

def align(query, target):
    '''
    Performs semi-global alignment using a variation of the Needleman-Wunsch
    algorithm and computes the associated cost matrix.
    :param query: read string
    :param target: genome string (at certain index and range)
    :return: cost matrix for alignment.
    '''

    M = np.zeros((len(query)+1, len(target)+1))

    for i in range(0, M.shape[0]):
        M[i][0] = i * gap

    for i in range(1, M.shape[0]):
        for j in range(1, M.shape[1]):
            cost = match
            if (query[i - 1] != target[j - 1]):
                cost = mismatch
            M[i][j] = min(M[i-1][j] + gap, M[i][j-1] + gap, M[i-1][j-1] + cost)

    return M

def backtrack(M, query, sequence):
    '''
    Performs backtracking given alignment matrix computed by align()
    to generate optimal alignment strings for given query and target string.
    :param M: Cost matrix for alignment returned by align()
    :param query: read string
    :param sequence: genome string (at certain index and range)
    :return: aligned strings
    '''

    endIndexI = M.shape[0] - 1
    endIndexJ = np.argmin(M[endIndexI])
    queryI = ''
    queryJ = ''
    #queryI = query[endIndexI - 1]
    #queryJ = sequence[endIndexJ - 1]
    while (endIndexI >= 1 and endIndexJ >= 1):

        if (M[endIndexI][endIndexJ] == M[endIndexI - 1][endIndexJ] + gap):
            queryI = query[endIndexI - 1] + queryI
            queryJ = '-' + queryJ
            endIndexI -= 1
        elif (M[endIndexI][endIndexJ] == M[endIndexI][endIndexJ - 1] + gap):
            queryI = '-' + queryI
            queryJ = sequence[endIndexJ - 1] + queryJ
            endIndexJ -= 1
        else:
            queryI = query[endIndexI - 1] + queryI
            queryJ = sequence[endIndexJ - 1] + queryJ
            endIndexI -= 1
            endIndexJ -= 1

    return queryI, queryJ, endIndexJ # - 1 # -1 because we added extra column/row for M


def cigar(query, target):
    '''
    Given aligned sequences, generate the CIGAR string.
    Supports matches =, mismatches X, insertion I, deletion D and
    soft clipping S.
    :param M: Cost matrix for alignment returned by align()
    :param query: aligned read string
    :param target: aligned genome string (at certain index and range)
    :return: CIGAR string for the aligned sequences
    '''

    cigar = ''

    match_count = 0
    mismatch_count = 0
    insertion_count = 0
    deletion_count = 0

    type = -1 # 0 mismatch, 1 match, 2 insertion, 3 deletion, 4 soft clipping

    posStart = 0
    posEnd = len(query) - 1

    while ((posStart < len(query)) & (query[posStart] == '-')):
        posStart += 1

    while ((posEnd >= 0) & (query[posEnd] == '-')):
        posEnd -= 1

    if posStart > 0:
        cigar += str(posStart)
        cigar += 'S'

    for i in range(posStart, posEnd + 1):

        if ((query[i] != target[i]) & (query[i] != '-') & (target[i] != '-')): # mismatch
            if type == 0:
                mismatch_count += 1
            elif type != 0:
                mismatch_count = 1
                type = 0
                if(match_count > 0):
                    cigar += str(match_count)
                    cigar += '='
                elif(insertion_count > 0):
                    cigar += str(insertion_count)
                    cigar += 'I'
                elif(deletion_count > 0):
                    cigar += str(deletion_count)
                    cigar += 'D'
                match_count = 0
                insertion_count = 0
                deletion_count = 0

        elif ((query[i] == target[i]) & (query[i] != '-') & (target[i] != '-')): # match
            if type == 1:
                match_count += 1
            elif type != 1:
                match_count = 1
                type = 1
                if(mismatch_count > 0):
                    cigar += str(mismatch_count)
                    cigar += 'X'
                elif(insertion_count > 0):
                    cigar += str(insertion_count)
                    cigar += 'I'
                elif(deletion_count > 0):
                    cigar += str(deletion_count)
                    cigar += 'D'
                mismatch_count = 0
                insertion_count = 0
                deletion_count = 0

        elif ((query[i] != target[i]) & (query[i] != '-') & (target[i] == '-')): # insertion
            if type == 2:
                insertion_count += 1
            elif type != 2:
                insertion_count = 1
                type = 2
                if(match_count > 0):
                    cigar += str(match_count)
                    cigar += '='
                elif(mismatch_count > 0):
                    cigar += str(mismatch_count)
                    cigar += 'X'
                elif(deletion_count > 0):
                    cigar += str(deletion_count)
                    cigar += 'D'
                match_count = 0
                mismatch_count = 0
                deletion_count = 0

        elif ((query[i] != target[i]) & (query[i] == '-') & (target[i] != '-')): # deletion
            if type == 3:
                deletion_count += 1
            elif type != 3:
                deletion_count = 1
                type = 3
                if(match_count > 0):
                    cigar += str(match_count)
                    cigar += '='
                elif(mismatch_count > 0):
                    cigar += str(mismatch_count)
                    cigar += 'X'
                elif(insertion_count > 0):
                    cigar += str(insertion_count)
                    cigar += 'I'
                match_count = 0
                mismatch_count = 0
                insertion_count = 0

    if(match_count > 0):
        cigar += str(match_count)
        cigar += '='

    elif(mismatch_count > 0):
        cigar += str(mismatch_count)
        cigar += 'X'

    elif(insertion_count > 0):
        cigar += str(insertion_count)
        cigar += 'I'

    elif(deletion_count > 0):
        cigar += str(deletion_count)
        cigar += 'D'

    if posEnd < (len(query) - 1):
        cigar += str(len(query) - posEnd - 1)
        cigar += 'S'

    return cigar

def getCigarFromIndex(query, target):
    M = align(query, target)
    x, y, index = backtrack(M, query, target)
    return cigar(x, y), np.min(M[M.shape[0] - 1]), index


'''
q = 'ATTCGAA'
t = 'ARUIAATCGAATTCA'
cigar, cost, index = getCigarFromIndex(q, t)
print(cigar)
print(cost)
print(index)
#print(np.min(M[M.shape[0] -1]))
#print(cigar(x, y))
#x = '-AT--TA--'
#y = 'AATCG-AGT'
#print(cigar(x,y))
#print(reverse_complement("ATCGTAAGT"))
'''
