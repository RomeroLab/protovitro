import numpy
from Bio.SubsMat import MatrixInfo


M = MatrixInfo.blosum62
Mrev = dict(((k[1],k[0]),M[k]) for k in M) # reverse the order of the AAs
M.update(Mrev)


def calc_bitscore(M,seq1,seq2,gap_opening=-11,gap_extension=-1):

    # the first position gets a little complicated
    # because for gap opening/extension we ask if pos-1 was a gap
    # solution: add a new character in front of sequences
    #  and only iterate from second postion
    seq1 = 'X'+seq1
    seq2 = 'X'+seq2

    score = 0
    for pos in range(1,len(seq1)): # from second position

        # both positions are not gaps
        if not(seq1[pos]=='-' and seq2[pos]=='-'): 

            # seq1 has a gap
            if seq1[pos]=='-':
                # seq1 previous postion was a gap: extension
                if seq1[pos-1]=='-':
                    score += gap_extension
                else: # opening and extension
                    score += gap_opening + gap_extension

            # seq2 has a gap
            elif seq2[pos]=='-':
                if seq2[pos-1]=='-':
                    score += gap_extension
                else:
                    score += gap_opening + gap_extension

            # seq1 and seq2 have a match
            else:
                score += M[(seq1[pos],seq2[pos])]

    # Karlin-Altschul parameters for ungapped alignments 
    # only apply for 11,1 gap penalties
    lam = 0.267
    logK = -3.194183212277829 # K = 0.041
    log2 = 0.69314718055994529

    bitscore = (lam*score - logK)/log2

    return bitscore
