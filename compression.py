# -*- coding: utf-8 -*-
# @Time    : 23.10.20
# @Authors  : Chenxi-Nie / Guillaume-Thiry
# @E-mails    ：chenie@student.ethz.ch / gthiry@student.ethz.ch
# @FileName: compression.py

## IMPORT

from input import read_fa
from bwt import Burrows_Wheeler_Transform

## FUNCTIONS

def save_string(string, filename):
    '''
    Save a string in a .txt file
    :param string: string to save
    :param filename: name of the save file (not the path)
    :return: Nothing
    '''

    directory = "A:/Mes documents (HDD)/CBM/Project 1/save/" #to change for each computer
    file = directory + filename + ".txt"
    f = open(file, 'w')
    f.write(string)
    f.close()

def rle_compression(seq):
    '''
    Does a basic Run-Lenght Encoding
    :param seq: sequence to compress (string)
    :return: the compressed sequence (string)
    '''

    n = len(seq)
    i = 1
    res = ''
    previous = seq[0]
    count = 1
    while(i<n):
        if seq[i] != previous:
            if count == 1:
                res += previous
            else:
                res += (str(count) + previous)
                count = 1
            previous = seq[i]
        else:
            count += 1
        i += 1
    if count == 1:
        res += previous
    else:
        res += (str(count) + previous)
        count = 1
    return res

def rle_decompression(seq):
    '''
    Does the Run-Lenght Encoding decompression
    :param seq: sequence to decompress (string)
    :return: the decompressed sequence (string)
    '''

    res = ""
    n = 0
    count = []
    for i in range(len(seq)):
        try:    #the character is a number, so we add that to count
            count.append(int(seq[i]))
        except: #the character is a letter, and we append the result
            if count != []:
                for k in range(len(count)): #convert the figures of the count into the real number
                    n *= 10
                    n += count[k]
                res += seq[i]*n
                count = []
                n = 0
            else:
                res += seq[i]
    return res


#Dicts used for the binary (de)compression

code_dict = {"A" : "000",
             "C" : "001",
             "G" : "010",
             "N" : "011",
             "T" : "100",
             "$" : "101"}

decode_dict = {"000" : "A",
             "001" : "C",
             "010" : "G",
             "011" : "N",
             "100" : "T",
             "101" : "$"}

def binary(s):
    '''
    Given a binary sequence, computes the corresponding decimal number
    :param s: binary sequence (string of 0s and 1s)
    :return: the decimal number (int)
    '''
    c = 1
    r = 0
    for i in range(len(s)-1,-1,-1):
        if s[i] == "1":
            r += c
        c *= 2
    return r

def debinary(num):
    '''
    Given a integer <127, computes its binary sequence (7 bits)
    :param num: integer to get into binary
    :return: the binary sequence (string)
    '''
    n = num
    k = 64
    res = ''
    while(k>=1):
        if n>=k:
            res += '1'
            n = n-k
        else:
            res += '0'
        k /= 2
    return res

def binary_compression(seq):
    '''
    Compresses a genome sequence using binary sequences and ASCII characters
    :param seq: genome sequence to compress (string)
    :return: compressed sequence (string)
    '''
    n = len(seq)
    i = 0
    current = ""
    res = ""
    while(i<n):
        char = seq[i]
        current += code_dict[char]
        if len(current)>=7:
            code = binary(current[0:7])
            current = current[7:]
            res += chr(code)
        i+=1
    res = current + '\n' + res
    return res

def binary_decompression(seq):
    '''
    Decompresses a genome sequence using binary sequences and ASCII characters
    :param seq: genome sequence to decompress (string)
    :return: decompressed sequence (string)
    '''
    n = len(seq)
    i = 0
    current = ""
    res = ""
    end = ""
    while(seq[i]!="\n"):
        end += seq[i]
        i += 1
    i += 1
    while(i<n):
        current += debinary(ord(seq[i]))
        i += 1
        while(len(current)>2):
            res += decode_dict[current[0:3]]
            current = current[3:]
    current = current + end
    while(len(current)>0):
        res += decode_dict[current[0:3]]
        current = current[3:]
    return res


## EXAMPLES
'''
fa_example = "A:/Mes documents (HDD)/CBM/Project 1/data_small/-CSCO-3h--genome.chr22.5K.fa"
genome = read_fa(fa_example)

bwt, rotated_genome = Burrows_Wheeler_Transform(genome)

rle = rle_compression(bwt)
bin = binary_compression(bwt)

print()
print("Original length : ", len(bwt))
print("RLE length : ", len(rle))
print("Binary length : ", len(bin))

print()
print(rle_decompression(rle) == bwt)
print(binary_decompression(bin) == bwt)
'''
