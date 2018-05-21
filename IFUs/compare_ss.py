#! /misc/local/python-2.7.11/bin/python

import sys, cPickle, string
from Bio import pairwise2
from parse_sequence_and_ss import PDBInfo

ss_dict = {'HT':'^',
           'HE':'*',
           'H ':'%',
           'E ':'%',
           'EH':'*',
           'ET':'*',
           ' T':'%',
           ' H':'%',
           ' E':'%',
           ' -':'-',
           'H-':'-',
           'E-':'-'}

j_scores = {'*':1.0,'%':0.6,'^':0.25,'-':0.0,'.':0.0}

minHammDiff = 0.5
minGraftSize = 10
minHammDist = minGraftSize*minHammDiff

def ss_string(s1,s2):

    s = ''

    for i in xrange(len(s1)):

        if s1[i] == s2[i]:
            s += '.'

        else:
            s += ss_dict[s1[i]+s2[i]]

    return s

def special_len(s):

    return(len([x for x in s if x != '-']))

def juicy_differences(ss_diff):

    diffs = string.split(ss_diff,'.')

    scores = []
    lengths = []
    spaces = []

    for i in xrange(len(diffs)):
        
        if i == 0:
            spaces.append(0)
        
        if not diffs[i]:
            spaces[-1] += 1
            continue

        scores.append(sum([j_scores[x] for x in diffs[i]]))
        lengths.append(special_len(diffs[i]))
        spaces.append(1)

    if diffs[-1]:
        spaces[-1] = 0
    else: 
        spaces[-1] -= 1

    return diffs,scores, lengths, spaces

def find_idx(iname,l):

    for i in xrange(len(l)):
        if l[i] == iname:
            return i 

def find_max_score(scoresN,lengthsN,spacesN,diffsN,name,scoresC,lengthsC,
                   spacesC,diffsC,max_string):

    #max_score = max(scores)
    #max_idx = find_idx(max_score,scores)

    new_score1 = -50
    new_score2 = -50
    size1 = 500
    size2 = 500
    cur_score = sum([j_scores[x] for x in max_string])


    diffsN = [x for x in diffsN if x]
    diffsC = [x for x in diffsC if x]

    #max_string = diffs[max_idx]

    #print scores
    #print diffs[max_idx]
    #print spaces
    #print lengths
    #sys.exit()

    o1 = 0
    o2 = len(diffsN)-1

    string1 = ''
    string2 = ''

    #print 'Hello1, ',max_string

    #print 'SpacesN, DiffsN: ',len(spacesN), len(diffsN)

    while (o1 < len(scoresC) or o2 >= 0):

        if o2>= 0:
            new_score1 = cur_score+scoresN[o2]
            size1 = special_len(max_string)+lengthsN[o2]+spacesN[o2+1]

            #print diffsN[o2],spacesN[o2+1]

            string1 = diffsN[o2]+'.'*spacesN[o2+1]+max_string

        else:
            new_score1 = -50
            size1 = 0
            string1 = '.'

        if o1 < len(scoresC):
            new_score2 = cur_score+scoresC[o1]
            size2 = special_len(max_string)+lengthsC[o1]+spacesC[o1]
            string2 = max_string+'.'*spacesC[o1]+diffsC[o1]

        else:
            new_score2 = -50
            size2 = 0
            string2 = '.'

        if new_score1 < minHammDiff*special_len(string1) and \
           new_score2 < minHammDiff*special_len(string2):
            return max_string, cur_score

        if new_score1/special_len(string1) >= new_score2/special_len(string2):
            cur_score = new_score1
            max_string = diffsN[o2]+'.'*spacesN[o2+1]+max_string
            o2 -= 1

        else:
            cur_score = new_score2
            max_string = max_string+'.'*spacesC[o1]+diffsC[o1]
            o1+= 1

        #print name, max_string, cur_score, diffs

    return max_string, cur_score

def find_seq_idxs(frag,seq):

    for i in xrange(len(seq)-len(frag)+1):
        if seq[i:i+len(frag)] == frag:
            return[i,i+len(frag)]

    return [-999,-999]

def get_discrepancies(ss_str):

    important_indices = []
    important_strings = []

    final_indices = []
    final_strings = []

    for i in xrange(len(ss_str)-minGraftSize):
        if sum([j_scores[x] for x in \
                ss_str[i:i+minGraftSize]])/float(minGraftSize) >= 0.5:
            important_indices.append(i)
            important_strings.append(ss_str[i:i+minGraftSize])
            
    elongating = 0
    cur_string = ''
    cur_X1 = -999
    cur_X2 = -999

    #print important_indices
    #print important_strings

    for i in xrange(len(important_indices)-1):

        #print important_indices[i]
        #print important_strings[i]

        if important_indices[i+1]-important_indices[i] == 1:
            #print 'Here1'
            if not elongating:
                elongating = 1
                cur_string = important_strings[i]+important_strings[i+1][-1]
                cur_X1 = important_indices[i]
                cur_X2 = important_indices[i+1]+minGraftSize

            else:
                cur_string = cur_string+important_strings[i+1][-1]
                cur_X2 = important_indices[i+1]+minGraftSize

        else:
            #print 'Here2'
            if elongating:
                final_indices.append((cur_X1,cur_X2))
                final_strings.append(cur_string)
            else:
                cur_string = important_strings[i]
                cur_X1 = important_indices[i]
                cur_X2 = important_indices[i]+minGraftSize

            elongating = 0

    if cur_X1 != -999 and cur_X2 != -999:

        final_indices.append((cur_X1,cur_X2))

    #print final_indices
    #print final_strings
            

    return final_indices


if __name__ == '__main__':

    ss_pred = open(sys.argv[1])
    pross   = open(sys.argv[2])

    SS = cPickle.load(ss_pred)
    PROSS = cPickle.load(pross)

    i = 0

    for k in SS.keys():

        #if k != '3KUYA':
        #    continue

        ss_seq = SS[k].sq
        ss_ss  = SS[k].ss

        try:

            ps_seq = PROSS[k].sq
            ps_ss  = PROSS[k].ss
            
        except KeyError:
            continue

        if ss_seq == '' or ps_seq == '':
            continue

        if ss_seq != ps_seq:
            #print 'Sequences not identical ', k

            alignment = pairwise2.align.localxs(ss_seq,ps_seq,-5000,-1000)
            #print alignment[0]
            continue

        else:
            ss_str = ss_string(ss_ss,ps_ss)

            discrepancies = get_discrepancies(ss_str)
            #print discrepancies
            
            for j in xrange(len(discrepancies)):
                
                #print ss_str


                #print ss_str[0:discrepancies[j][0]]
                
                diffsN,scoresN, lengthsN, spacesN = juicy_differences(ss_str[0:discrepancies[j][0]])


                #print diffsN,scoresN, lengthsN, spacesN

                #print ss_str[discrepancies[j][1]:]

                diffsC,scoresC,lengthsC,spacesC = juicy_differences(ss_str[discrepancies[j][1]:])

                #print  diffsC,scoresC,lengthsC,spacesC

                #print 'Hello, ',ss_str[discrepancies[j][0]:discrepancies[j][1]]

                max_string, max_score = find_max_score(scoresN,lengthsN,
                                                       spacesN,diffsN,k,
                                                       scoresC,lengthsC,
                                                       spacesC,diffsC,
                                                       ss_str[discrepancies[j][0]:discrepancies[j][1]])

                if special_len(max_string) > minGraftSize:
               
                    idxs = find_seq_idxs(max_string,ss_str)
                    
                    #if not ss_seq[idxs[0]:idxs[1]]:
                        

                    print k, max_string, max_score, max_score/special_len(max_string),ss_seq[idxs[0]:idxs[1]]
                    #print [x for x in diffsN if x],diffsN,scoresN, lengthsN, spacesN
                    #print [x for x in diffsC if x],diffsC,scoresC,lengthsC,spacesC
                    #print ss_str
                    #print 'Hello, ',ss_str[discrepancies[j][0]:discrepancies[j][1]]
                    #print idxs[0],idxs[1], len(ss_seq)

                    #print ss_seq
                    #print ss_str
            

                #print ss_str[idxs[0]:idxs[1]],max_string
                #print ss_ss[idxs[0]:idxs[1]]
                #print ps_ss[idxs[0]:idxs[1]]
                #print ss_seq[idxs[0]:idxs[1]]
                #print ps_seq[idxs[0]:idxs[1]]

        i += 1

        #if i == 5:
        #    sys.exit()
