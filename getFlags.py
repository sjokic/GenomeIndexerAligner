# Authors: Stefan Jokic, Adrian Taubner

def getFlags(Read1_isFwd, Read2_isFwd, isMapped1, isMapped2):
    flag1 = 1 + 64
    flag2 = 1 + 128

    if(isMapped1 and isMapped2):
        flag1 += 2
        flag2 += 2

    elif(isMapped1 and (not isMapped2)):
        flag1 += 8
        flag2 += 4

    elif((not isMapped1) and isMapped2):
        flag1 += 4
        flag2 += 8

    elif((not isMapped1) and (not isMapped2)):
        flag1 += 4 + 8
        flag2 += 4 + 8

    if(Read1_isFwd and (not Read2_isFwd)):
        flag1 += 32
        flag2 += 16
    elif((not Read1_isFwd) and Read2_isFwd):
        flag1 += 16
        flag2 += 32

    return flag1, flag2
