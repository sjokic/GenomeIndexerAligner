# -*- coding: utf-8 -*-
# @Time    : 29.10.20 11:21
# @Author  : Guillaume T.
# @FileName: bwt_linear_time.py

# Implementing a Suffix Array construction in linear time
# Reference article : Linear Work Suffix Array Construction
# By : Karkkainen, Sanders, Burkhardt


# Comparison of pairs
def comp_pairs(a,b):
    if (a[0]<b[0]):
        return True
    elif ((a[0]==b[0]) and (a[1]<=b[1])):
        return True
    else:
        return False

# Comparison of triplets
def comp_triplets(a,b):
    if (a[0]<b[0]):
        return True
    elif (a[0]==b[0]):
        return comp_pairs((a[1],a[2]), (b[1],b[2]))
    else:
        return False

# One radix-sort pass (offset indicates which digit is used)
def radix(ind, table, n, K, offset=0):

    c = [0]*(K+1)   #Alphabet counter

    for i in range(n): # Counting the characters
        c[table[ind[i]+offset]] += 1

    positions = []
    sum = 0
    for i in range(K+1):
        positions.append(sum)
        sum += c[i]

    res = [0]*len(ind)

    for i in range(n):
        res[positions[table[ind[i]+offset]]] = ind[i]
        positions[table[ind[i]+offset]] += 1

    return res

def getI(SA12, t, n0):
    if SA12[t]<n0:
        return SA12[t]*3+1
    else:
        return (SA12[t]-n0)*3+2

def Suffix_Array(T, n, K):

    n0 = int((n+2)/3)
    n1 = int((n+1)/3)
    n2 = int(n/3)
    n02 = n0+n2

    #Step 0
    R = []
    for i in range(n+(n0-n1)):
        if i%3 != 0:
            R.append(i)
    R += [0,0,0]


    #Step 1

    R1 = radix(R,T,n02,K,2)
    R2 = radix(R1,T,n02,K,1)
    R3 = radix(R2,T,n02,K,0)


    current_rank = 0
    current = "///"

    for i in range(n02):
        triplet = str(T[R3[i]]) + str(T[R3[i]+1]) + str(T[R3[i]+2])
        if triplet != current:
            current = triplet
            current_rank += 1
        if R3[i]%3 == 1:
            R[int(R3[i]/3)] = current_rank
        else:
            R[int(R3[i]/3) + n0] = current_rank


    if current_rank < n02:
        SA12 = Suffix_Array(R, n02, current_rank)
        for i in range(n02):
            R[SA12[i]] = i+1
    else:
        SA12 = [0]*(n02+3)
        for i in range(n02):
            SA12[R[i]-1] = i
    #Step 2

    R0 = []
    for i in range(n02):
        if SA12[i]<n0:
            R0.append(3*SA12[i])
    SA0 = radix(R0, T, n0, K)

    #Step 3
    t = n0-n1
    p = 0
    k = 0
    SA = []
    while(k<n):
        i = getI(SA12, t, n0)
        j = SA0[p]

        if SA12[t]<n0:
            comp = comp_pairs((T[i],R[SA12[t]+n0]), (T[j], R[int(j/3)]))
        else:
            comp = comp_triplets((T[i],T[i+1],R[SA12[t]-n0+1]),(T[j],T[j+1],R[int(j/3)+n0]))


        if comp:
            SA.append(i)
            t += 1
            if(t==n02):
                k+=1
                while(p<n0):
                    SA.append(SA0[p])
                    p+=1
                    k+=1
        else:
            SA.append(j)
            p+=1
            if(p==n0):
                k+=1
                while(t<n02):
                    SA.append(getI(SA12,t,n0))
                    t+=1
                    k+=1
        k+=1
    return SA
