#! /misc/local/python-2.7.11/bin/python

import sys, cPickle

from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.SubsMat import MatrixInfo as matlist

class PDBInfo:

    def __init__(self,block):

        self.ID = block[0]
        self.sq = block[1]
        self.ss = block[2]

similarity_dict = {'A':{'A':1,'C':0,'D':0,'E':0,'F':0,'G':0,'H':0,'I':0.5,'K':0,'L':0.5,'M':0.5,'N':0,'P':0,
                        'Q':0,'R':0,'S':0,'T':0,'V':0.5,'W':0,'Y':0},
                   'C':{'A':0,'C':1,'D':0,'E':0,'F':0,'G':0,'H':0,'I':0,'K':0,'L':0,'M':0,'N':0,'P':0,
                        'Q':0,'R':0,'S':0.5,'T':0,'V':0,'W':0,'Y':0},
                   'D':{'A':0,'C':0,'D':1,'E':0.75,'F':0,'G':0,'H':0.5,'I':0,'K':0.5,'L':0,'M':0,'N':0.5,'P':0,
                        'Q':0.5,'R':0.5,'S':0.5,'T':0.5,'V':0,'W':0,'Y':0},
                   'E':{'A':0,'C':0,'D':0.75,'E':1,'F':0,'G':0,'H':0.5,'I':0,'K':0.5,'L':0,'M':0,'N':0.5,'P':0,
                        'Q':0.5,'R':0.5,'S':0.5,'T':0.5,'V':0,'W':0,'Y':0},
                   'F':{'A':0,'C':0,'D':0,'E':0,'F':1,'G':0,'H':0,'I':0.5,'K':0,'L':0.5,'M':0,'N':0,'P':0,
                        'Q':0,'R':0,'S':0,'T':0,'V':0,'W':0,'Y':0.5},
                   'G':{'A':-0.5,'C':-0.5,'D':-0.5,'E':-0.5,'F':-0.5,'G':1,'H':-0.5,'I':-0.5,'K':-0.5,'L':-0.5,
                        'M':-0.5,'N':-0.5,'P':-0.5,'Q':-0.5,'R':-0.5,'S':-0.5,'T':-0.5,'V':-0.5,'W':-0.5,
                        'Y':-0.5},
                   'H':{'A':0,'C':0,'D':0.5,'E':0.5,'F':0,'G':0,'H':1,'I':0,'K':0.5,'L':0,'M':0,'N':0.5,'P':0,
                        'Q':0.5,'R':0.5,'S':0.5,'T':0.5,'V':0,'W':0,'Y':0},
                   'I':{'A':0.5,'C':0,'D':0,'E':0,'F':0.5,'G':0,'H':0,'I':1,'K':0,'L':0.75,'M':0.5,'N':0,'P':0,
                        'Q':0,'R':0,'S':0,'T':0,'V':0.75,'W':0,'Y':0},
                   'K':{'A':0,'C':0,'D':0.5,'E':0.5,'F':0,'G':0,'H':0.5,'I':0,'K':1,'L':0,'M':0,'N':0.5,'P':0,
                        'Q':0.5,'R':0.75,'S':0.5,'T':0.5,'V':0,'W':0,'Y':0},
                   'L':{'A':0.5,'C':0,'D':0,'E':0,'F':0.5,'G':0,'H':0,'I':0.75,'K':0,'L':1,'M':0.5,'N':0,'P':0,
                        'Q':0,'R':0,'S':0,'T':0,'V':0.75,'W':0,'Y':0},
                   'M':{'A':0.5,'C':0,'D':0,'E':0,'F':0,'G':0,'H':0,'I':0.5,'K':0,'L':0.5,'M':1,'N':0,'P':0,
                        'Q':0,'R':0,'S':0,'T':0,'V':0.5,'W':0,'Y':0},
                   'N':{'A':0,'C':0,'D':0.5,'E':0.5,'F':0,'G':0,'H':0.5,'I':0,'K':0.5,'L':0,'M':0,'N':1,'P':0,
                        'Q':0.75,'R':0.5,'S':0.5,'T':0.5,'V':0,'W':0,'Y':0},
                   'P':{'A':0.5,'C':0,'D':0,'E':0,'F':0,'G':0,'H':0,'I':0.5,'K':0,'L':0,'M':0,'N':0,'P':1,
                        'Q':0,'R':0,'S':0,'T':0,'V':0,'W':0,'Y':0},
                   'Q':{'A':0,'C':0,'D':0.5,'E':0.5,'F':0,'G':0,'H':0.5,'I':0,'K':0.5,'L':0,'M':0,'N':0.75,'P':0,
                        'Q':1,'R':0.5,'S':0.5,'T':0.5,'V':0,'W':0,'Y':0},
                   'R':{'A':0,'C':0,'D':0.5,'E':0.5,'F':0,'G':0,'H':0.5,'I':0,'K':0.75,'L':0,'M':0,'N':0.5,'P':0,
                        'Q':0.5,'R':1,'S':0.5,'T':0.5,'V':0,'W':0,'Y':0},
                   'S':{'A':0,'C':0,'D':0.5,'E':0.5,'F':0,'G':0,'H':0.5,'I':0,'K':0.5,'L':0,'M':0,'N':0.5,'P':0,
                        'Q':0.5,'R':0.5,'S':1,'T':0.75,'V':0,'W':0,'Y':0},
                   'T':{'A':0,'C':0,'D':0.5,'E':0.5,'F':0,'G':0,'H':0.5,'I':0,'K':0.5,'L':0,'M':0,'N':0.5,'P':0,
                        'Q':0.5,'R':0.5,'S':0.75,'T':1,'V':0,'W':0,'Y':0},
                   'V':{'A':0.5,'C':0,'D':0,'E':0,'F':0.5,'G':0,'H':0,'I':0.75,'K':0,'L':0.75,'M':0.5,'N':0,
                        'P':0,'Q':0,'R':0,'S':0,'T':0,'V':1,'W':0,'Y':0},
                   'W':{'A':0,'C':0,'D':0,'E':0,'F':0.5,'G':0,'H':0,'I':0.25,'K':0,'L':0.25,'M':0,'N':0,'P':0,
                        'Q':0,'R':0,'S':0,'T':0,'V':0,'W':1,'Y':0},
                   'Y':{'A':0,'C':0,'D':0,'E':0,'F':0.5,'G':0,'H':0,'I':0.25,'K':0,'L':0.25,'M':0,'N':0,'P':0,
                        'Q':0,'R':0,'S':0,'T':0,'V':0,'W':0,'Y':1}}


def get_pdb_info_block(f,idx):

    pdbID = f[idx][1:5]+f[idx][6]

    #print pdbID

    j = idx
    ss_found = 0
    
    seq = ''
    ss  = ''

    while j < len(f) and (f[j][0] != '>' or f[j][1:5] + f[j][6] == pdbID):
        if not ss_found and f[j][0] != '>':
            seq += f[j]
        if ss_found:
            ss += f[j]
        if not ss_found and f[j][0] == '>' and f[j][-6:] == 'secstr':
            ss_found = 1
        
        j += 1

        if j < len(f) and not f[j]:
            j += 1
    

    return j, [pdbID, seq, ss]

def string_similarity(s1,s2):

    return sum([s1[i]==s2[i] for i in xrange(len(s1))])

def sequence_equivalency(s1,s2):
    

    return sum([similarity_dict[s1[i]][s2[i]] for i in xrange(len(s2))])

def good_fold_switch(seq1,seq2,ss1,ss2):

    short_seq = ''
    long_seq  = ''
    short_ss  = ''
    long_ss   = ''

    if len(seq1) > len(seq2):
        short_seq = seq2
        short_ss = ss2
        long_seq = seq1
        long_ss = ss1
    else:
        short_seq = seq1
        short_ss  = ss1
        long_seq = seq2
        long_ss  = ss2

    for i in xrange(len(long_seq)-len(short_seq)+1):
        
        ss_similarity = string_similarity(short_ss,long_ss[i:i+len(short_seq)])/float(len(short_seq))

        if ss_similarity <= 0.25:
            
            equivalency_score=sequence_equivalency(long_seq[i:i+len(short_seq)],short_seq)
            if equivalency_score/float(len(short_seq)) >= 0.4:

                print equivalency_score
                print short_seq
                print short_ss
                print long_ss[i:i+len(short_seq)]
                print long_seq[i:i+len(short_seq)]


if __name__ == '__main__':

    f= open(sys.argv[1]).read().splitlines()
    
    of = open(sys.argv[2],'w')

    i = 0
    
    pdb_dict = {}

    while i < len(f):

        i, block = get_pdb_info_block(f,i)

        #print i, block

        pdb_dict[block[0]] = PDBInfo(block)
    
    cPickle.dump(pdb_dict,of)

#good_fold_switch(pdb_dict['1OSPO'].sq, pdb_dict['2V64E'].sq[7:33],pdb_dict['1OSPO'].ss,pdb_dict['2V64E'].ss[7:33])












