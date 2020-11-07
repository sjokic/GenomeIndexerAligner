# Author: Adrian Taubner

import os

from input import read_sam


def compareResult(place, result1File, result2File): # result1File is correct file
    result1path = place + result1File
    result2path = place + result2File

    print(result1path)
    result1pos, result1cigar = read_sam(result1path, True)
    result2pos, result2cigar = read_sam(result2path, False)

    totalFalse = 0

    for i in range(0, len(result1pos)):
        if result1cigar[i] == result2cigar[i]:
            a=0
            #print("index: " + str(i) + ", correct")
        else:
            print("index: " + str(i) + ", incorrect, correct Pos: " + str(result1pos[i]) + ", our Pos:" + str(result2pos[i]) +
                  ", correct cigar: " + str(result1cigar[i]) + ", our cigar: " + str(result2cigar[i]))
            totalFalse += 1

    print("Total false = " + str(totalFalse))

result1 = "output_tiny_30xCov.sam"
result2 = "our_output_tiny_30xCov.sam"
folder = os.getcwd() + "/data_small/"

compareResult(folder, result1, result2)
