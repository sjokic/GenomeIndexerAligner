# Authors: Stefan Jokic, Adrian Taubner

from align import reverse_complement

def writeToFile(f, header, sn, ln, reads1, reads2):
    if(header):
        versionNumber = 1.4
        s0 = 'unsorted'
        # Header
        f.write('@HD' + '\t')
        f.write('VN:' + str(versionNumber) + '\t')
        f.write('SO:' + s0 + '\n')
        f.write('@SQ' + '\t')
        f.write('SN:' + sn + '\t')
        f.write('LN:' + str(ln) + '\n')

    mapq = 0
    rnext = '='
    pnext = 0
    tlen = 0
    seq = '*'
    qual = '*'

    # alignment section
    seq1 = ''
    seq2 = ''
    qual1 = ''
    qual2 = ''
    if(reads1[5]):
        seq1 = reads1[4]
        qual1 = reads1[6]
    else:
        seq1 = reverse_complement(reads1[4])
        qual1 = reads1[6][::-1]

    if(reads2[5]):
        seq2 = reads2[4]
        qual2 = reads2[6]
    else:
        seq2 = reverse_complement(reads2[4])
        qual2 = reads2[6][::-1]

    f.write(reads1[0] + '\t' + str(reads1[1]) + '\t' + sn + '\t' + str(reads1[2]) + '\t' + str(mapq) + '\t' + reads1[3] + '\t'
        + rnext + '\t' + str(reads2[2]) + '\t' + str(tlen) + '\t' + seq1 + '\t' + qual1 + '\n')
    f.write(reads2[0] + '\t' + str(reads2[1]) + '\t' + sn + '\t' + str(reads2[2]) + '\t' + str(mapq) + '\t' + reads2[3] + '\t'
        + rnext + '\t' + str(reads1[2]) + '\t' + str(tlen) + '\t' + seq2 + '\t' + qual2 + '\n')
