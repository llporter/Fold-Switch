#! /misc/local/python-2.7.11/bin/python

import sys, string
import Bio.PDB
from Bio import pairwise2
import urllib
from difflib import SequenceMatcher
from fuzzywuzzy import fuzz
from LINUS import evaluator, molecule
from biomol import pdbparse
from biomol.pdbout import pack_pdb_line

PDBDIR = ''
MIN_RMSD = 4.0
STR_SIM = 0.8
MIN_LEN = 12

aa3to1 = {'ALA':'A','CYS':'C','ASP':'D','GLU':'E','PHE':'F','GLY':'G','HIS':'H',\
          'ILE':'I','LYS':'K','LEU':'L','MET':'M','ASN':'N','PRO':'P','GLN':'Q',\
          'ARG':'R','SER':'S','THR':'T','VAL':'V','TRP':'W','TYR':'Y','MSE':'X',\
          'NAG':'','BMA':'','MAN':'','SO4':'','NHE':'', ' CL':'', 'HOH':'',
          ' CU':'', ' ZN':'','BNG':'','PLM':'','ACT':'','ATM':'',' MG':'', 'UNK':'X'}#,'TYD':''}

AAs = ['ALA','CYS','ASP','GLU','PHE','GLY','HIS','ILE','LYS','LEU','MET','ASN',\
       'PRO','GLN','ARG','SER','THR','VAL','TRP','TYR']

class PDBINFO:

    def __init__(self,info,indices,sequence,ss):

        self.ID          = None
        self.chain       = None
        self.length      = None
        self.description = None
        self.indices     = None
        self.sequence    = None
        self.ss          = None

        self.get_pdb_info(info,indices,sequence,ss)

        self.resis2include = None

    def get_pdb_info(self,info,indices,sequence,ss):

        self.ID = info[0][1:5]
        self.chain = info[0][-1]

        #print 'info ', len(info), info

        if len(info) == 1:
            
            self.length = None
            self.description = None

        else:
            self.length = int(info[2][7:])
            self.description = string.join(info[3:],' ')
            
        self.indices = (int(indices.split()[1]),int(indices.split()[3]))
        self.sequence = sequence
        self.ss = ss

        #self._print()

    def _print(self):

        print 'ID= ',self.ID
        print 'Chain= ',self.chain
        print 'length= ',self.length
        print 'description= ',self.description
        print 'indices= ', self.indices
        print 'sequence= ', self.sequence
        print 'ss= ', self.ss

def max_str(diffstring,maxDiff,sequence):

    l = len(maxDiff)

    if diffstring==maxDiff:
        return sequence

    for i in xrange(len(diffstring)-l):

        if maxDiff == diffstring[i:i+l]:
            return sequence[i:i+l]

def min_str(str1,str2):

    s1 = string.join([x for x in str1 if x!='-'],'');
    s2 = string.join([x for x in str2 if x!='-'],'');

    if len(s1)<=len(s2):
        return s1

    return s2

class LOG3:

    def __init__(self,f):

        self.template   = None
        self.pdbInfo    = []
        self.identities = []
        self.hammDists  = []
        self.hammDiffs  = []
        self.maxSeqs    = []
        self.maxStrings = []
        self.diffStrings= []
        self.error1 = 0

        self.parse(f)

        self.bseqs = []
        
        self.fss = []
        self.rfs = []
        self.fps = []

    def parse(self,fname):
        
        i = 0

        f = open(fname).read().splitlines()

        while i < len(f):

            info = f[i].split()
            #print info

            if info[0] == 'Template:':
                self.template = PDBINFO(f[i+2].split()[:-1],f[i+5],f[i+3],f[i+4])
                i += 8
                #self.template._print()
                continue
            
            if info[0][0] == '>':
                
                pdb1 = PDBINFO(info,f[i+3],f[i+1],f[i+2])

                sequence = ''

                #pdb1._print()
                
                self.pdbInfo.append(pdb1)
                self.identities.append(SequenceMatcher(None,self.template.sequence,
                                                       pdb1.sequence).ratio())
                
                self.hammDists.append(float(f[i+4].split()[2][:-1]))
                self.hammDiffs.append(float(f[i+4].split()[5][:-1]))
                self.maxStrings.append(f[i+4].split()[-4][:-1])
                self.diffStrings.append(f[i+5].split()[0])

                sequence = max_str(self.diffStrings[-1],self.maxStrings[-1],pdb1.sequence)

                self.maxSeqs.append((f[i+4].split()[-1],sequence))
                
                i+=8

            else:
                print 'File parsing error %s' %(fname)
                self.error1 = 1
                i += 1

        #print self.maxSeqs

    def pdb2idx(self,info):

        for i in xrange(len(info)):

            if info[i] == '&' and info[i+1][0]=='>':
                return i

    def get_pdb(self,ID,chain):
        
        pdb_parser = Bio.PDB.PDBParser(QUIET = True)

        try:
            pdb = pdb_parser.get_structure(ID,PDBDIR+ID+'.pdb')
            
        except IOError:
            fetch_pdb(ID)
            pdb = pdb_parser.get_structure(ID,PDBDIR+ID+'.pdb')

        try:
            return pdb[0][chain]
        except:
            return 0

    def remove_breaks(self,seq,ss):

        fseq = ''

        for i in xrange(len(ss)):
            if ss[i] != '-':
                fseq = fseq + seq[i]

        return fseq

    def revised_indices(self,pdb,seq,ss):

        sidx1 = None
        eidx1 = None

        keys = pdb.child_dict.keys()
        
        resis = []

        for k in keys:
            resis.append(k[1])

        resis.sort()

        wseq = string.join([aa3to1[x.get_resname()] for x in pdb],'')

        seq = self.remove_breaks(seq,ss)

        #print seq
        
        l = len(seq)

        if len(wseq)==l and seq == wseq:
            sidx1 = resis[0]
            eidx1 = resis[0]+l

        elif len(wseq) < l:
            sidx1 = resis[0]
            eidx1 = resis[0]+len(wseq)

        else:
            for i in xrange(len(wseq)-l):
                if seq ==  wseq[i:i+l]:
                    sidx1 = resis[0]+i
                    eidx1 = resis[0]+i+l

        #print sidx1, eidx1

        return [sidx1,eidx1]

        
        
    def get_segment(self,pdb,indices,seq,ss,good_residues):

        ref = []
            
        if not good_residues:
           # print 'Here'
            good_residues = range(indices[0],indices[1]+1)

        try:
            ref.append(pdb[good_residues[0]]['CA'])

        except KeyError:
            #print "here2"

            keys = pdb.child_dict.keys()
        
            resis = []

            for k in keys:
                resis.append(k[1])

            resis.sort()
            
            good_residues = [x+resis[0] for x in good_residues]

            #print good_residues
            #sys.exit()

        for i in good_residues:
            #print pdb[i].get_resname()
            ref.append(pdb[i]['CA'])
                  
        return ref

    def redundant(self,idx,descriptions):

        redundant = 0

        for d in descriptions:

            if fuzz.partial_ratio(self.pdbInfo[idx].description,d) >= STR_SIM:
                redundant = 1

            if redundant:
                return 1

        return redundant

    def get_resi_offset(self,pdb):

        keys = pdb.child_dict.keys()
        
        resis = []

        for k in keys:
            resis.append(k[1])

        resis.sort()

        return resis[0]


    def correct_indices(self,idx,pdb1,pdb2):

        usequence1 = None
        usequence2 = None

        long_seq = None
        short_seq = None

        good_seq_frags1 = None
        good_seq_frags1 = None

        good_indices1 = []
        good_indices2 = []

        #print self.maxSeqs[idx]

        sequence1 = self.maxSeqs[idx][0]

        if self.template.sequence==None or sequence1 == None:
            return -998,-998

        ss1       = max_str(self.template.sequence,sequence1,self.template.ss)
        indices1  = self.template.indices

        sequence2 = self.maxSeqs[idx][1]

        if self.pdbInfo[idx].sequence == None or sequence2 == None:
            return -998,-998

        ss2       = max_str(self.pdbInfo[idx].sequence,sequence2,self.pdbInfo[idx].ss)
        indices2  = self.pdbInfo[idx].indices

        #if '-' not in ss2 and '-' not in ss1:
            #return 0 ,0

        #print self.maxSeqs[idx]

        usequence1 = self.remove_breaks(sequence1,ss1)
        usequence2 = self.remove_breaks(sequence2,ss2)

        #print usequence1, usequence2

        if len(usequence1) < len(usequence2):
            short_seq = usequence1
            long_seq  = usequence2
        else:
            short_seq = usequence2
            long_seq  = usequence1

        alignments = pairwise2.align.localxs(short_seq,long_seq,-5000,-1000)

        good_seq_frags = string.split(alignments[0][0],'-')
        good_seq_frags = [x for x in good_seq_frags if x != '']

        #print good_seq_frags

        good_seq_frags2 = string.split(alignments[0][0],'-')
        good_seq_frags2 = [x for x in good_seq_frags2 if x != '']

        #wseq1 = string.join([aa3to1[x.get_resname()] for x in pdb1],'')
        #wseq2 = string.join([aa3to1[x.get_resname()] for x in pdb2],'')

        #print good_seq_frags, good_seq_frags2

        wseq1 = get_sequence(pdb1)
        wseq2 = get_sequence(pdb2)

        wa = pairwise2.align.localxs(wseq1,wseq2,-5000,-1000)

        #print wa

        wa1=wa[0][0][wa[0][-2]:wa[0][-1]]
        wa2=wa[0][1][wa[0][-2]:wa[0][-1]]

        for i in xrange(len(good_seq_frags)):

            if len(good_seq_frags[i]) < 4:
                continue
            
            a = pairwise2.align.localxs(good_seq_frags[i],wa1,-5000,-1000)

            minseq1 = min_str(a[0][0],a[0][1])

            a = pairwise2.align.localxs(minseq1,wseq1,-5000,-1000)

            ab = string.join([a[0][0][x] for x in xrange(len(a[0][0]))\
                              if a[0][0][x] != '-' and a[0][1][x] != '-' ],'')

            #print "Here: ",ab,wseq1

            if ((ab not in wseq1) and fuzz.partial_ratio(ab,wseq1) < 70) \
               or len(ab) < MIN_LEN:
                a = -999
                a2 = -999
                break
            
            #good_indices1 += [x for x in xrange(len(a[0][0])) if a[0][0][x] != '-' and \
                              #a[0][1][x] != '-']

            a2 = pairwise2.align.localxs(good_seq_frags[i],wa2,-5000,-1000)

            minseq2 = min_str(a2[0][0],a2[0][1])

            a2  = pairwise2.align.localxs(minseq2,wseq2,-5000,-1000)

            #print a2

            ab2 = string.join([a2[0][0][x] for x in xrange(len(a2[0][0]))\
                              if a2[0][0][x] != '-' and a2[0][1][x] != '-' ],'')

            #ab2 = string.join([x for x in a2[0][0] if x != '-'],'')

            #print "Here: ",ab2,wseq2,fuzz.partial_ratio(ab2,wseq2)

            if ((ab2 not in wseq2) or fuzz.partial_ratio(ab2,wseq2) <70) \
               or len(ab2) < MIN_LEN:
                a2 = -999
                break

            #good_indices2 += [x for x in xrange(len(a2[0][0])) if a2[0][0][x] != '-' and \
                              #a2[0][1][x] != '-']

            if len(ab2) != len(ab):

                final_alignment = pairwise2.align.localxs(ab,ab2,-5000,-1000)

                if final_alignment[0][2] < MIN_LEN:
                    a = -999
                    a2 = -999
                    break

                final_seq = string.join([final_alignment[0][0][x] for x in xrange(len(final_alignment[0][0]))\
                                         if final_alignment[0][0][x] != '-' and final_alignment[0][1][x] != '-'],'')
                
                a1 = pairwise2.align.localxs(final_seq,wseq1,-5000,-1000)
                good_indices1 += [x for x in xrange(len(a1[0][0])) if a1[0][0][x] != '-' and \
                                  a1[0][1][x] != '-']

                a2 = pairwise2.align.localxs(final_seq,wseq2,-5000,-1000)
                good_indices2 += [x for x in xrange(len(a2[0][0])) if a2[0][0][x] != '-' and \
                                  a2[0][1][x] != '-']

            else:
               good_indices1 += [x for x in xrange(len(a[0][0])) if a[0][0][x] != '-' and \
                                 a[0][1][x] != '-'] 
               good_indices2 += [x for x in xrange(len(a2[0][0])) if a2[0][0][x] != '-' and \
                                 a2[0][1][x] != '-']
                
            
        #print string.join([wseq1[x] for x in good_indices1])
        #print string.join([wseq2[x] for x in good_indices2])

       # print wseq2

        #print a1, a2



        if a2 == -999 or a == -999:
            return a,a2

        #return [x+indices1[0] for x in good_indices1], [x+indices2[0] for x in good_indices2]

        #offset1 = self.get_resi_offset(pdb1)
        #offset2 = self.get_resi_offset(pdb2)
        
        #good_indices1 = [x+offset1 for x in good_indices1]
        #good_indices2 = [x+offset2 for x in good_indices2]

        #good_indices1 = [pdb1[x].id for x in good_indices1]
        #good_indices2 = [pdb2[x].id for x in good_indices2]

        ids1 = get_ids(pdb1)
        ids2 = get_ids(pdb2)

        good_ids1 = []
        good_ids2 = []

        for i in xrange(len(ids1)):
            if i in good_indices1:
                good_ids1.append(ids1[i])

        for i in xrange(len(ids2)):
            if i in good_indices2:
                good_ids2.append(ids2[i])
                
        

        #return good_indices1, good_indices2

        return good_ids1, good_ids2
        
    def verify(self):

        pdb1 = ''
        pdb2 = ''

        switch_descriptions = []
        fp_descriptions = []

        switch_found = 0
        false_positive = 0

        bad_indices = open("bad_indices.txt",'a')
        residue_insertions = open("bad_insertions.txt",'a')
        

        super_imposer = Bio.PDB.Superimposer()

        for i in xrange(len(self.pdbInfo)):

            #print self.template.ID+self.template.chain,\
                  #self.pdbInfo[i].ID+self.pdbInfo[i].chain

            #If redundant hit, document and continue

            if len(switch_descriptions) > 0 and self.redundant(i,switch_descriptions):
                #print '%s and %s are redundant fold switches.  Skipping...' \
                    #%(self.template.ID +self.template.chain,
                      #self.pdbInfo[i].ID +self.pdbInfo[i].chain)
                self.rfs.append((self.template.ID +self.template.chain,
                                 self.pdbInfo[i].ID +self.pdbInfo[i].chain))
                continue

            #If redundant false positive, document and continue

            if len(fp_descriptions) > 0 and self.redundant(i,fp_descriptions):
                #print '%s and %s are redundant false positives.  Skipping...' \
                    #%(self.template.ID+self.template.chain,
                      #self.pdbInfo[i].ID+self.pdbInfo[i].chain)
                self.fps.append((self.template.ID+self.template.chain,
                                 self.pdbInfo[i].ID+self.pdbInfo[i].chain))
                continue

            pdb1 = self.get_pdb(self.template.ID, self.template.chain)

            if not pdb1:
                print 'Invalid pdb %s%s.  Skipping...' %(self.template.ID,
                                                         self.template.chain)
                continue

            insertions1 = find_insertions(pdb1,self.template.ID+self.template.chain)
            
            pdb2 = self.get_pdb(self.pdbInfo[i].ID, self.pdbInfo[i].chain)

            if not pdb2:
                print 'Invalid pdb %s%s.  Skipping...' %(self.pdbInfo[i].ID,
                                                         self.pdbInfo[i].chain)
                continue

            insertions2 = find_insertions(pdb2,self.pdbInfo[i].ID+self.pdbInfo[i].chain)
            
            indices1,indices2 = self.correct_indices(i,pdb1,pdb2)

            self.template.resis2include = indices1
            self.pdbInfo[i].resis2include = indices2

            if indices1 == -999:
                
                print "Missing electron density in %s%s.  Skipping..." %(self.template.ID,
                                                                         self.template.chain)
                continue

            if indices2 == -999:

                print "Missing electron density in %s%s.  Skipping..." %(self.pdbInfo[i].ID,
                                                                         self.pdbInfo[i].chain)
                continue

            if indices1 == -998:

                print 'Alignment too short between %s%s and %s%s.  Skipping...' %(self.template.ID,
                                                                                  self.template.chain,
                                                                                  self.pdbInfo[i].ID,
                                                                                  self.pdbInfo[i].chain)
                continue

            if len(indices1) != len(indices2):
                bad_indices.write("Bad Indices: %s %s\n"\
                                  %(self.template.ID+self.template.chain,\
                                    self.pdbInfo[i].ID+self.pdbInfo[i].chain))
                continue
                
            
            pdb_ref1 = self.get_segment(pdb1,self.template.indices,
                                        self.template.sequence,self.template.ss,
                                        self.template.resis2include)
            
            pdb_ref2 = self.get_segment(pdb2,self.pdbInfo[i].indices,
                                        self.pdbInfo[i].sequence,self.pdbInfo[i].ss,
                                        self.pdbInfo[i].resis2include)

            super_imposer.set_atoms(pdb_ref1,pdb_ref2)
            super_imposer.apply(pdb2.get_atoms())

            #print super_imposer.rms

            if super_imposer.rms >= MIN_RMSD:
                #print 'Length: ', len(indices1), 'RMSD: ',super_imposer.rms
                #print self.template.ID+self.template.chain,\
                    #self.pdbInfo[i].ID+self.pdbInfo[i].chain
                switch_found = 1
                switch_descriptions.append(self.template.description)
                switch_descriptions.append(self.pdbInfo[i].description)

                #print switch_descriptions
                    
                #print 'Switch found from %s to %s!' %(self.template.ID +\
                                                      #self.template.chain,
                                                      #self.pdbInfo[i].ID+\
                                                      #self.pdbInfo[i].chain)

                #seq1 = string.join([aa3to1[i.get_resname()] for i in pdb_ref1],'')
                #seq2 = string.join([aa3to1[i.get_resname()] for i in pdb_ref2],'')

                self.fss.append((self.template.ID +self.template.chain,
                                 self.pdbInfo[i].ID+self.pdbInfo[i].chain,
                                 super_imposer.rms))

                try:

                    Bio.PDB.Dice.extract(pdb1,self.template.chain,
                                         self.template.resis2include[0][1],
                                         self.template.resis2include[-1][1],
                                         'frags/'+self.template.ID+self.template.chain+\
                                         '_'+str(self.template.resis2include[0][1])+'_'+\
                                         str(self.template.resis2include[-1][1])+'.pdb')

                
                    Bio.PDB.Dice.extract(pdb2,self.pdbInfo[i].chain,
                                         self.pdbInfo[i].resis2include[0][1],
                                         self.pdbInfo[i].resis2include[-1][1],
                                         'frags/'+self.pdbInfo[i].ID+self.pdbInfo[i].chain+\
                                         '_'+str(self.pdbInfo[i].resis2include[0][1])+'_'+\
                                         str(self.pdbInfo[i].resis2include[-1][1])+'.pdb')

                except TypeError:
                    residue_insertions.write("Bad Insertions.  No structures: %s %s\n"\
                                      %(self.template.ID+self.template.chain,\
                                        self.pdbInfo[i].ID+self.pdbInfo[i].chain))
                    
                    
            else:
                fp_descriptions.append(self.template.description)
                fp_descriptions.append(self.pdbInfo[i].description)
                
                #print '%s and %s are a false_positive.'  %(self.template.ID +\
                                                           #self.template.chain,
                                                           #self.pdbInfo[i].ID+\
                                                           #self.pdbInfo[i].chain)
                self.fps.append((self.template.ID +self.template.chain,
                                 self.pdbInfo[i].ID+self.pdbInfo[i].chain))

        bad_indices.close()
        residue_insertions.close()
                
                
        return 1

    def write_info(self,id,list,f):

        for l in list:
            f.write('%10s %10s %10s\n' %(id,l[0],l[1]))
            
    def write_info2(self,id,list,f):

        for l in list:
            f.write('%10s %10s %10s %10.4f\n' %(id,l[0],l[1],l[2]))
            
    def out(self,id,fps,rfs,fss):

        self.write_info(id,self.fps,fps)
        self.write_info(id,self.rfs,rfs)
        self.write_info2(id,self.fss,fss)

    def out2(self,id):

        for l in self.fss:
            print 'Fold switch: %10s %10s %10s %10.4f' %(id,l[0],l[1],l[2])

        
def fetch_pdb(id):
    url = 'http://www.rcsb.org/pdb/files/%s.pdb' % id
    return urllib.urlretrieve(url,PDBDIR+'%s.pdb' %id)

def remove_hetero(chain):

    for res in chain:
        name = res.get_resname()
        id = res.id
        if name not in AAs:
            chain.detach_child(id)

    return chain

def get_sequence(chain):

    s = ''

    for res in chain:
        if res.id[0] == ' ':
            s += aa3to1[res.get_resname()]

    return s

def find_insertions(chain,pdbID):

    insertion_ids = []

    for res in chain:
        if res.id[2] != ' ':
            insertion_ids.append(res.id)
            print 'Insertion: ',pdbID, res.id

    return insertion_ids

def get_ids(chain):

    ids = []

    for res in chain:
        if res.id[0] == ' ':
            ids.append(res.id)

    return ids
            
                

if __name__ == '__main__':

    #f = open(sys.argv[1]).read().splitlines()
    
    f = sys.argv[1]

    #n= string.split(string.split(sys.argv[2],'.')[0],'_')[1]
    #n= "test"

    n = sys.argv[2]
    #fps = open("false_positives_"+n+"_2.txt",'a')
    #rfs = open("redundant_fold_switches_"+n+"_2.txt",'a')
    #fss = open("unique_fold_switches_"+n+"_2.txt",'a')

    #for i in xrange(len(f)):

    fold_switch_pairs = LOG3(f)

    if fold_switch_pairs.error1:
        sys.exit()

    real_switches = fold_switch_pairs.verify()

    fold_switch_pairs.out2(f)

    #fps.close()
    #rfs.close()
    #fss.close()

        
        
        
