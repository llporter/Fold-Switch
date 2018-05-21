#! /groups/looger/home/porterl/python/bin/python

import sys, cPickle, string
from frag_mvalu import *
from verify_log3s_coop import *

###Change Variables below to match locations on your workstation.
###Also make sure to change the python path above to match where your version of python is installed.
###NOTE: THIS SOFTWARE RUNS ON PYTHON 2.7; other versions of python will give errors.

#Location of PDB from which Independent Folding Units will be calculated
pdbPath = ''
#Extension of PDB File
pdbExtension = '.pdb'
#If you have a directory that contains the whole PDB, insert here
PDBDIR = '/groups/looger/loogerlab/DSSP/PDBs/'
#Where you want your output files to be saved
FILEDIR3 = '/groups/looger/home/porterl/m_value_calculations/fs_110917/'

#Minimum Qualifying Ratio (QR) for a fold switch.  0.78 was the value used in Porter ahd Looger, PNAS 2018.
MIN_RAT = 0.78

aa3to1 = {'ALA':'A','CYS':'C','ASP':'D','GLU':'E','PHE':'F','GLY':'G','HIS':'H',\
          'ILE':'I','LYS':'K','LEU':'L','MET':'M','ASN':'N','PRO':'P','GLN':'Q',\
          'ARG':'R','SER':'S','THR':'T','VAL':'V','TRP':'W','TYR':'Y','MSE':'X',\
          'NAG':'X','BMA':'X','MAN':'X'}

class PDBInfo:
    
    def __init__(self,block):
        
        self.ID = block[0]
        self.sq = block[1]
        self.ss = block[2]

def generate_clean_pdb(pdbID,chain):

    fetch_pdb(pdbID)

    chn = get_protein_chain(pdbID,chain,pdbPath,pdbExtension)
    try:
        clean_chain(chn)
    except:
        'Protein chain %s not found.  Skipping...' %(fnames[i])
        return 0
    
    writePDB(chn,pdbID+chain)
    return 1

def idx4ID2(FSPs,ID2):

    for i in xrange(len(FSPs.pdbInfo)):

        if FSPs.pdbInfo[i].ID+FSPs.pdbInfo[i].chain == ID2:
            return i
    
    return -999

def reconstruct_string(SS_strings,space_ends,SS):

    s = ''

    for i in xrange(len(SS_strings)):

        s += SS_strings[i]
        l = len(s)
        if i < len(space_ends):
            if SS_strings[0] != 0:
                s += ' '*(space_ends[i]-l)
            else:
                s += ' '*(space_ends[i+1]-l)


def ss_boundaries(SS):

    SS_strings = SS.split()

    space_ends = []

    s = -999

    SS_map = []

    for i in xrange(len(SS)):
        if SS[i] == ' ':
            s = i
            
        elif i == 0 or SS[i-1] ==' ':
            s = i

        SS_map.append(s)

    return SS_map

def ss2pos(SS_strings,ss):

    finished = 0
    i = 0
    j = 0

    pos = []
    cpos = []
    good_ss = []

    while not finished:
        while i < len(ss):

            if j == len(SS_strings):
                finished = 1
                break

            
            if 'H' not in SS_strings[j] and\
               'E' not in SS_strings[j]:
                j+= 1
                i+= 1
                continue

            if SS_strings[j] == ss[i:i+len(SS_strings[j])]:
                good_ss.append(SS_strings[j])
                pos.append(i)
                i += len(SS_strings[j])
                cpos.append(i)
                j += 1
                continue

            i += 1

    posDict = {}

    for i in xrange(len(pos)):
        posDict[pos[i]] = i

    return pos, cpos, posDict
            
            

def smart_map(ss):

    SS_strings = SS.split()

    s = -999

    sMap = []

    pos, posDict = ss2pos(SS_strings,ss)
        

def getfSeqLoc(seq,fSeq):
    
    alignment = pairwise2.align.localxs(fSeq,seq,-5000,-1000)#

    if alignment[0][1][0] == '-':
        
        o = 0
        for i in xrange(len(alignment[0][1])):
            if alignment[0][1][i] == '-':
                o += 1
            if alignment[0][1][i] != '-':
                break

        print alignment[0]

        print seq[alignment[0][3]-o:alignment[0][4]-o]

        return [alignment[0][3]-o,alignment[0][4]-o]
            
        
    return [alignment[0][3],alignment[0][4]]

def ss_boundariesC(SS_mapN):

    j = 0
    sm1 = ''

    mapC = []

    for i in xrange(len(SS_mapN)):
        
        if SS_mapN[i] != sm1:
            sm1 = SS_mapN[i]
            mapC += [i-1]*(i-j)
            j = i

    l = len([x for x in SS_mapN if x == SS_mapN[-1]])

    if l == 1:
        mapC.append(SS_mapN[-1])
    else:
        mapC += [SS_mapN[-1]+l-1]*l

    return mapC

def getMol(name):

    saveout = sys.stdout
    sys.stdout = open(os.devnull, 'w')
    m = molecule.MoleculeWithBackboneHydrogens(name)
    sys.stdout = saveout

    return m

def get_segment_Mvals(pdbID,segment,pdbin):

    cosolvents = ['Urea']

    resi_dict = get_resi_breaks(pdbin)
    
    ACE = []
    nt_offset=1
    if resi_dict[segment[0]] == 0:
        print 'Warning: could not add ACE blocking group.  M-value '+\
            'calculation will begin at second residue in segment.'
    else:
        for atm in pdbin[resi_dict[segment[0]-1]+1:\
                         resi_dict[segment[0]-1]+4]:
            if atm[12:16] == ' CA ':
                atm = atm[:12]+' CH3'+atm[16:]
                
            ACE.append(atm[:17]+'ACE '+atm[21:])

    NME = []

    if segment[1] >= len(resi_dict.keys()):
        segment[1] = len(resi_dict.keys())-1
    if resi_dict[segment[1]] == len(pdbin):
        ct_offset = 0
        print 'Warning: could not add NME blocking group.  M-value '+\
            'calculation will end at second to last residue in segment.'
    else:
        for atm in pdbin[resi_dict[segment[1]]:\
                         resi_dict[segment[1]]+2]:
            if atm[12:16] == ' CA ':
                atm = atm[:12]+' CH3'+atm[16:]

            NME.append(atm[:17]+'NME '+atm[21:])
            
    segpdb =ACE+pdbin[resi_dict[segment[0]]:resi_dict[segment[1]]]+NME

    
    segmentMvalues = get_total_mvalues(segpdb,cosolvents,[1,1])
        
    return segmentMvalues


def printMvalInfo(cosolvents,segmentMvalues,pdbMvalues):

    for c in cosolvents:
                segSum = sum(segmentMvalues.mvalues[c])
                situSum= sum(pdbMvalues.mvalues[c])
                print segSum, situSum

                BBDiff = -sum(pdbMvalues.folded_BB_ASAs)+\
                             sum(segmentMvalues.folded_BB_ASAs)
                SCDiff = -sum(pdbMvalues.folded_SC_ASAs)+\
                             sum(segmentMvalues.folded_SC_ASAs)


                mval_ratio = segSum/situSum

                bbsegMvals = sum(segmentMvalues.bbmvals[c])
                scsegMvals = sum(segmentMvalues.scmvals[c])
                bbpdbMvals = sum(pdbMvalues.bbmvals[c])
                scpdbMvals = sum(pdbMvalues.scmvals[c])

                for i in xrange(len(segmentMvalues.bbmvals[c])):
                    print segmentMvalues.bbmvals[c][i], pdbMvalues.bbmvals[c][i], segmentMvalues.scmvals[c][i], pdbMvalues.scmvals[c][i]


                print sum(pdbMvalues.folded_BB_ASAs),sum(segmentMvalues.folded_BB_ASAs), sum(pdbMvalues.folded_SC_ASAs), sum(segmentMvalues.folded_SC_ASAs)
                print ((bbsegMvals+scsegMvals)/(bbpdbMvals+scpdbMvals))
                print bbsegMvals, bbpdbMvals, scsegMvals, scpdbMvals
                print bbpdbMvals+scpdbMvals
                

def calculate_ratio2(cosolvents,segmentMvalues,pdbMvalues):

    for c in cosolvents:
        segSum = sum(segmentMvalues.mvalues[c])
        situSum= sum(pdbMvalues.mvalues[c])

    return segSum/situSum

def calculate_ratio(cosolvents,segmentMvalues,pdbMvalues,segment):

    for c in cosolvents:
        segSum = sum(segmentMvalues.mvalues[c])
        situSum= sum(pdbMvalues.mvalues[c][segment[0]:segment[1]])

    return segSum/situSum

def removeSeqDashes(seq):

    return string.join([x for x in seq if x not in ['-',' ']],'')

def closestn(lst,n):

    mindiff = 999
    minidx = -1

    for i in xrange(len(lst)):

        print n, lst[i]

        if i==0 and n-lst[i] < 0:
            return n

        if n-lst[i] > 0 and n-lst[i] < mindiff:
            mindiff = n-lst[i]
            minidx = i

        if i >0 and n-lst[i] <0:
            break

    return lst[minidx]

def closestc(lst,n):

    mindiff = 999
    minidx = -1

    for i in xrange(len(lst)):

        print n, lst[i]

        if lst[i]-n > 0 and lst[i]-n < mindiff:
            mindiff = lst[i]-n
            minidx = i

        if lst[i]-n > 0:
            break

    return lst[minidx]

def expand_secondary_structure(pdbID,segment,pdbin,pdbMvalues,SS_mapN,SS_mapC,nres,maxSeqLocs,
                               molSeq,sequence,cosolvents,ssPos,sscPos,posDict,rnumDict):

    ratios = []

    nterm = None
    cterm = None

    print sequence[maxSeqLocs[0]:maxSeqLocs[1]]

    flen = len(sequence[maxSeqLocs[0]:maxSeqLocs[1]])

    #used to be segment[0]
    if posDict.has_key(maxSeqLocs[0]):
        nterm = maxSeqLocs[0]
    else:
        nterm = closestn(ssPos,maxSeqLocs[0])
        print ssPos, nterm

    #used to be segment[1]
    if posDict.has_key(maxSeqLocs[1]):
        cterm = maxSeqLocs[1]
    else:
        cterm = closestc(ssPos,maxSeqLocs[1])
        print ssPos, cterm

    print sequence[nterm:cterm]
    print molSeq[rnumDict[nterm]:rnumDict[cterm]]


    iran = range(max(0,posDict[nterm]-3),posDict[nterm]+1)
    jran = range(posDict[cterm],min(len(ssPos),posDict[cterm]+4))


    leni = len(iran)
    lenj = len(jran)

        

    for i in xrange(max(0,posDict[nterm]-3),posDict[nterm]+1):
        for j in xrange(posDict[cterm],min(len(ssPos),posDict[cterm]+4)):

            if ssPos[i] > ssPos[j]:
                continue
        
            segment = [rnumDict[ssPos[i]],rnumDict[ssPos[j]]]


            segmentMvalues = get_segment_Mvals(pdbID,segment,pdbin)

            r = calculate_ratio(cosolvents,segmentMvalues,pdbMvalues,[segment[0],min(nres,segment[1]+1)])

            ratios.append(r)

    print ratios

    maxi = -999
    maxj = -999
    maxk = -999

    maxr = -999

    good_ranges = []
            
    for k in xrange(len(ratios)):
        if ratios[k] >= 0.7:
            good_ranges.append([ratios[k],rnumDict[ssPos[iran[k/lenj]]], rnumDict[ssPos[jran[k%lenj]]]])

    if good_ranges:
        return good_ranges
    else:
        return[[-999,-999,-999]]


def expand_secondary_structure2(pdbID,segment,pdbin,pdbMvalues,SS_mapN,SS_mapC,nres,maxSeqLocs,
                               molSeq,sequence,cosolvents,ssPos,sscPos,posDict):


    ratios = []

    nterm = None
    cterm = None

    print sequence[maxSeqLocs[0]:maxSeqLocs[1]]

    if posDict.has_key(segment[0]):
        nterm = segment[0]
    else:
        nterm = closestn(ssPos,segment[0])
        print ssPos, nterm

    if posDict.has_key(segment[1]):
        cterm = segment[1]
    else:
        cterm = closestc(ssPos,segment[1])
        print ssPos, cterm


    print 'Expand test'
    print segment

    print sequence[maxSeqLocs[0]:maxSeqLocs[1]]
    print molSeq[nterm:cterm]

    if not posDict.has_key(nterm):
        posDict[nterm] = 0
        for i in xrange(len(ssPos)):
            posDict[ssPos[i]] += 1

        ssPos = [nterm]+ssPos

    print ssPos

    iran = range(max(0,posDict[nterm]-3),posDict[nterm]+1)
    jran = range(posDict[cterm],min(len(ssPos),posDict[cterm]+4))

    leni = len(iran)
    lenj = len(jran)

        

    for i in xrange(max(0,posDict[nterm]-3),posDict[nterm]+1):
        for j in xrange(posDict[cterm],min(len(ssPos),posDict[cterm]+4)):

            if ssPos[i] > ssPos[j]:
                continue
        

            segment = [ssPos[i],ssPos[j]]


            segmentMvalues = get_segment_Mvals(pdbID,segment,pdbin)

            r = calculate_ratio(cosolvents,segmentMvalues,pdbMvalues,[segment[0],min(nres,segment[1]+1)])

            ratios.append(r)

    print ratios

    maxi = -999
    maxj = -999
    maxk = -999

    maxr = -999

    good_ranges = []

    for k in xrange(len(ratios)):
        if ratios[k] >= 0.7:
            good_ranges.append([ratios[k],ssPos[iran[k/lenj]], ssPos[jran[k%lenj]]])

    if good_ranges:
        return good_ranges
    else:
        return[[-999,-999,-999]]

    
def removeDoubleDashes(ss,seq):

    ss2 = ''
    sq2 = ''

    for i in xrange(len(ss)):
        if ss[i] != '-':
            ss2 += ss[i]
            sq2 += seq[i]

    return ss2, sq2

def tweak_boundaries(pdbID, segment, pdbin, nres, pdbMvalues,cosolvents,molSeq,
                     nstep,o,rorig,sorig):

    for i in xrange(nstep):
        if segment[0]+o-i < 0:
            break
        for j in xrange(nstep):

            seg = [max(segment[0]+o-i,0),min(segment[1]-o+j,nres)]

            segmentMvalues = get_segment_Mvals(pdbID,seg,pdbin)

            r = calculate_ratio(cosolvents,segmentMvalues,pdbMvalues,seg)

            print r, i, j, molSeq[seg[0]:seg[1]]

            if r >= MIN_RAT and r <= 10.0:
                 return r, seg

            if segment[1]-o+j >= nres:
                break

        if segment[0]+o-i <= 0:
            break

    return rorig,sorig 

def maximize_boundaries(pdbID, segment, pdbin, nres, pdbMvalues,cosolvents,molSeq,
                     nstep,o,rorig,sorig):


    maxRat = rorig
    maxSeg = sorig

    for i in xrange(nstep):
        if segment[0]+o-i < 0:
            break
        for j in xrange(nstep):

            seg = [max(segment[0]+o-i,0),min(segment[1]-o+j,nres)]

            segmentMvalues = get_segment_Mvals(pdbID,seg,pdbin)

            r = calculate_ratio(cosolvents,segmentMvalues,pdbMvalues,seg)

            print r, i, j, molSeq[seg[0]:seg[1]]

            if r > maxRat and r < 10.0:
                maxRat = r
                maxSeg = seg

            if segment[1]-o+j >= nres:
                break

        if segment[0]+o-i <= 0:
            break

    return maxRat, maxSeg

def tweak_boundaries_plus(pdbID, segment, pdbin, nres, pdbMvalues,cosolvents,molSeq,
                          nstep,o,rorig,sorig):

    maxR = -999
    maxS = []

    for j in xrange(nstep):

        seg = [max(segment[0],0),min(segment[1]-o+j,nres)]

        segmentMvalues = get_segment_Mvals(pdbID,seg,pdbin)

        r = calculate_ratio(cosolvents,segmentMvalues,pdbMvalues,seg)

        print r, j, molSeq[seg[0]:seg[1]]

        if r >= maxR:
            maxR = r
            maxS = seg

        if r >= MIN_RAT and r <= 10.0:
            return r, seg

        if segment[1]-o+j >= nres:
            break


    if maxR > rorig:
        return maxR, maxS

    return rorig,sorig 

def tweak_boundaries_minus(pdbID, segment, pdbin, nres, pdbMvalues,cosolvents,molSeq,
                          nstep,o,rorig,sorig):

    maxR = -999
    maxS = []

    for j in xrange(nstep):

        seg = [max(segment[0]-o-j,0),min(segment[1],nres)]

        segmentMvalues = get_segment_Mvals(pdbID,seg,pdbin)

        r = calculate_ratio(cosolvents,segmentMvalues,pdbMvalues,seg)

        print 'Here ', r, j
        print r, j, molSeq[seg[0]:seg[1]]

        if r >= maxR:
            maxR = r
            maxS = seg

        if r >= MIN_RAT and r <= 10.0:
            return r, seg

        if segment[0]+o-j <= 0:
            break

    if maxR > rorig:
        return maxR, maxS

    return rorig,sorig 

def add_terminal_indices(rnd,nres):

    k = rnd.keys()
    k.sort()

    j = 0

    if k[0] != 0:
        for i in xrange(0,k[0]):
            rnd[i] = 0

    if k[-1] < nres-1:
        for i in xrange(k[-1],nres):
            rnd[i] = rnd[k[-1]]

    return rnd
        
def make_rnumDict(molSeq,refSeq):

    alignment = pairwise2.align.localxs(molSeq,refSeq,-5,-0.1)

    molA = alignment[0][0]
    refA = alignment[0][1]

    print molA
    print refA
    print len(refSeq)

    rnumDict = {}

    j = 0
    
    print alignment

    for i in xrange(len(molA)):

        if molA[i] != '-':
            rnumDict[i] = j
            j += 1

    k = rnumDict.keys()
    k.sort()
    curIdx = -999
    curIdxp1 = -999
    j = 0
        
    for i in xrange(len(k)):

        if i == 0:
            curIdx = k[i]
            continue
            
        if k[i] == curIdx+1 or (k[i] == curIdx+j and j != 0):
            curIdx = k[i]
            if i < len(k)-1:
                curIdxp1 = k[i+1]
            j = 0

        else:
            while curIdx + j < k[i]:
                j += 1
                rnumDict[curIdx+j] = rnumDict[curIdxp1]

            j += 1

    add_terminal_indices(rnumDict,len(refSeq))

    return rnumDict
            

    
    
def pdb_cooperativity(pdbSequence,pdbSS,ID,chain,maxSeq):

    cosolvents = ['Urea']

    print pdbSequence
    print pdbSS

    ss,sequence = removeDoubleDashes(pdbSS,pdbSequence)

    pos, cpos, posDict = ss2pos(ss.split(),ss)

    sequence = removeSeqDashes(sequence)
    mSeq = removeSeqDashes(maxSeq)

    print mSeq
    
    generated = generate_clean_pdb(ID, chain)
    
    if not generated: return 0

    SS_mapN = ss_boundaries(ss)

    SS_mapC = ss_boundariesC(SS_mapN)

    maxSeqLocs = getfSeqLoc(sequence,mSeq)

    print maxSeqLocs

    if maxSeqLocs[1] >= len(sequence):
        maxSeqLocs[1] = len(sequence)-1
 
    o1=1

    if SS_mapN[maxSeqLocs[0]]-1 <= 0:
        o1 = 0

    o2=1
    
    if SS_mapC[maxSeqLocs[1]] >= len(sequence)-1:
        o2 = 0

    maxSSLocs = [SS_mapN[maxSeqLocs[0]]-o1,SS_mapC[maxSeqLocs[1]]+o2]

    if maxSSLocs[0] <= 5:
        maxSSLocs[0] = 0

    if maxSSLocs[1] >= len(sequence)-5:
        maxSSLocs[1] = len(sequence)

    maxSq = sequence[maxSSLocs[0]:maxSSLocs[1]]
    
    maxSq = removeSeqDashes(maxSq)

    m = getMol(ID+chain+'.pdb')

    molSeq = string.join([aa3to1[x] for x in m.resnam],'')

    nres = m.numres

    del m

    mol_rnumDict = make_rnumDict(molSeq,sequence)

    segment = [mol_rnumDict[maxSeqLocs[0]],mol_rnumDict[maxSeqLocs[1]]]

    pdbin = open(ID+chain+'.pdb').read().splitlines()


    pdbMvalues = get_total_mvalues(pdbin,cosolvents,[0,0])

    r, seg = tweak_boundaries(ID+chain,[segment[0],segment[1]],
                              pdbin,nres,pdbMvalues,cosolvents,molSeq,50,0,0.0,segment)

    if r >= MIN_RAT:
        r,seg = maximize_boundaries(ID+chain,[seg[0],seg[1]],
                                    pdbin,nres,pdbMvalues,cosolvents,molSeq,15,0,r,seg)


    return r, seg[0],seg[1], molSeq[seg[0]:seg[1]+1]

    if r2 >= 0.7 and r < 0.87:
        r, seg = tweak_boundaries_plus(ID+chain,seg2,
                                       pdbin,nres,pdbMvalues,cosolvents,molSeq,50,0,r2,seg2)

    if r >=0.87 or r2 >=0.87:
        if r2 >= 0.87 and r < 0.87:
            r = r2
            seg = seg2

        elif r >= 0.87 and r2 < 0.87:
            r = r
            seg = seg

        elif seg2[1]-seg2[0] < seg[1]-seg[0]:
            r = r2
            seg = seg2

    else:
        seg = [seg2[0],seg[1]]
        segmentMvalues = get_segment_Mvals(ID+chain,seg,pdbin)
        r = calculate_ratio(cosolvents,segmentMvalues,pdbMvalues,seg)

    good_ranges = []
    
    if r < 0.7 and segment[1]-segment[0] >=15:

        good_ranges = []

        r,seg = tweak_boundaries(ID+chain,[segment[0],segment[1]],
                                 pdbin,nres,pdbMvalues,cosolvents,molSeq,10,0,
                                 r,seg)

        if r < MIN_RAT:
        
            good_ranges = expand_secondary_structure(ID+chain,segment,pdbin,pdbMvalues,
                                                     SS_mapN,SS_mapC,nres,maxSeqLocs,molSeq,sequence,
                                                     cosolvents,
                                                     pos,cpos, posDict,mol_rnumDict)

        if (good_ranges and good_ranges[0][0] != -999):

            for i in xrange(len(good_ranges)):
                if r < 0.9:
                    r,seg = tweak_boundaries(ID+chain,[good_ranges[i][1],good_ranges[i][2]],
                                             pdbin,nres,pdbMvalues,cosolvents,molSeq,10,0,r,seg)
                if r >= MIN_RAT:
                    segment = seg
                    break

        elif r >= 0.75 and not good_ranges:
            
            if r < 0.9:

                r,seg = tweak_boundaries(ID+chain,[seg[0],seg[1]],
                                         pdbin,nres,pdbMvalues,cosolvents,molSeq,10,0,r,seg)

            else:
                if seg:
                    segment = seg
            
            
    if r >= 0.7 and r < 0.9 and good_ranges and good_ranges[0][0] != -999:

        r,seg = tweak_boundaries(ID+chain, segment, pdbin,nres,pdbMvalues,cosolvents,molSeq,10,0,r,seg)

    if not seg:
        seg = segment

    return r, seg[0],seg[1], molSeq[seg[0]:seg[1]+1]
    

    
if __name__ == '__main__':

    f = open(sys.argv[1]).read().splitlines()

    of_name = FILEDIR3+string.split(sys.argv[1],'.')[0]+'.fsc'

    DSSP_PDBs = cPickle.load(open(sys.argv[2]))

    of = open(of_name,'w')

    for i in xrange(len(f)):

        info1 = f[i].split()

        ID = info1[0][:4]
        chn = info1[0][4]
        seq = info1[-1]


        info = DSSP_PDBs[ID+chn]

        p1_cooperativity = pdb_cooperativity(info.sq,info.ss,ID,chn,seq)

        print info1[0], p1_cooperativity

        if p1_cooperativity[0] >= MIN_RAT:

            of.write('%5s %1.4f %5i %5i %20s %12s\n' %(ID+chn,p1_cooperativity[0],
                                                    p1_cooperativity[1],
                                                    p1_cooperativity[2],
                                                    p1_cooperativity[3],
                                                    seq))
        
                       
