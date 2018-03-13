import gzip, sys, string
from biomol import pdbparse

coil_path = '/databases/coil/nonwind/';
tor_path  = '/databases/torsions/pdb';
pdb_path  = '/databases/pdb/pdb';

motif_list = ['T1_','T1p','T2_','T2p','T3_','T3p','HA_','HB1','HB2','HB3',
              'HB4','HB5','HB6','HC_','HD1','HD2','HD3','HD4','HE1','HE2',
              'HF_','HG1','HG2','HG3','HG4','HG5','HG6','HG7','HG8','HG9',
              'HGA','HGB','HGC','HGD','HGE','HGF','HGG','IGT','P2A','P2B',
              'AL_']

all_motifs = ['IGT','T1_','T1p','T2_','T2p','T3_','T3p','HA_','HB1','HB2',
              'HB3','HB4','HB5','HB6','HC_','HD1','HD2','HD3','HD4','HE1',
              'HE2','HF_','HG1','HG2','HG3','HG4','HG5','HG6','HG7','HG8',
              'HG9','HGA','HGB','HGC','HGD','HGE','HGF','HGG','P2A','P2B',
              'AL_','DPA','A2D','HLX','EXT','CIS','H__','CL_']

motifs = {'IGT':(((-100, -70),(  50, 100)),),
          'T1_':((( -90, -30),( -60,   0)), ((-120, -60),( -30,  30))),
          'T1p':(((  30,  90),(   0,  60)), ((  60, 120),( -30,  30))),
          'T2_':((( -90, -30),(  90, 150)), ((  50, 110),( -30,  30))),
          'T2p':(((  30,  90),(-150, -90)), ((-110, -50),( -30,  30))),
          'T3_':((( -90, -30),( -60,   0)), (( -90, -30),( -60,   0))),
          'T3p':(((  30,  90),(   0,  60)), ((  30,  90),(   0,  60))),
          'HA_':(((-114, -90),( 115, 165)), (( -98, -74),( 125, 175))),
          'HB1':(((-115, -85),(  -5,  25)), ((-116, -88),(  -5,  25)),
                 ((-127, -97),( 106, 154))),
          'HB2':(((-115, -85),(  -5,  25)), ((-116, -88),( 155, 180)),
                 ((-127, -97),( 106, 154))),
          'HB3':(((-115, -85),(  -5,  25)), ((-116, -88),(-180,-175)),
                 ((-127, -97),( 106, 154))),
          'HB4':(((-115, -85),( 115, 145)), ((-116, -88),(  -5,  25)),
                 ((-127, -97),( 106, 154))),
          'HB5':(((-115, -85),( 115, 145)), ((-116, -88),( 155, 180)),
                 ((-127, -97),( 106, 154))),
          'HB6':(((-115, -85),( 115, 145)), ((-116, -88),(-180,-175)),
                 ((-127, -97),( 106, 154))),
          'HC_':(((-110, -54),( -62, -38)), ((-104, -76),(  70, 100)),
                 ((-107, -83),( -14,  14)), ((-127, -97),(  75, 125))),
          'HD1':((( -95, -79),( -34,  -6)), ((  76,  94),(   5,  35)),
                 ((-101, -79),( -56,  -4))),
          'HD2':((( -95, -79),( -34,  -6)), ((  76,  94),(   5,  35)),
                 ((-101, -79),( 127, 161))),
          'HD3':((( -95, -79),( 151, 169)), ((  76,  94),(   5,  35)),
                 ((-101, -79),( -56,  -4))),
          'HD4':((( -95, -79),( 151, 169)), ((  76,  94),(   5,  35)),
                 ((-101, -79),( 127, 161))),
          'HE1':((( -98, -82),( -15,  15)), ((  53,  67),(  28,  52)),
                 ((-102, -78),(  88, 112)), ((-100, -70),( -42, -18)),
                 ((-112, -68),(  90, 150))),
          'HE2':((( -98, -82),( -15,  15)), ((  53,  67),(  28,  52)),
                 ((-102, -78),( -22,   2)), (( -99, -71),(-144,-116)),
                 ((-112, -68),(  90, 150))),
          'HF_':((( -99, -81),( -16,  -2)), ((  79,  91),(   3,  21)),
                 ((-113, -87),( 110, 144)), ((-112, -88),( 118, 142))),
          'HG1':(((-143,-117),(  73,  97)), (( -69, -51),( -29, -15)),
                 ((-104, -76),( -32, -12)), ((-113, -87),( -40, -14)),
                 ((-104, -76),( -44, -16))),
          'HG2':(((-143,-117),(  73,  97)), (( -69, -51),( -29, -15)),
                 ((-104, -76),( -32, -12)), ((-113, -87),( -40, -14)),
                 ((-104, -76),( 107, 133))),
          'HG3':(((-143,-117),(  73,  97)), (( -69, -51),( -29, -15)),
                 ((-104, -76),( -32, -12)), ((-113, -87),( 108, 136)),
                 ((-104, -76),( -44, -16))),
          'HG4':(((-143,-117),(  73,  97)), (( -69, -51),( -29, -15)),
                 ((-104, -76),( -32, -12)), ((-113, -87),( 108, 136)),
                 ((-104, -76),( 107, 133))),
          'HG5':(((-143,-117),(  73,  97)), (( -69, -51),( -29, -15)),
                 ((-104, -76),(  20,  40)), ((-113, -87),( -40, -14)),
                 ((-104, -76),( -44, -16))),
          'HG6':(((-143,-117),(  73,  97)), (( -69, -51),( -29, -15)),
                 ((-104, -76),(  20,  40)), ((-113, -87),( -40, -14)),
                 ((-104, -76),( 107, 133))),
          'HG7':(((-143,-117),(  73,  97)), (( -69, -51),( -29, -15)),
                 ((-104, -76),(  20,  40)), ((-113, -87),( 108, 136)),
                 ((-104, -76),( -44, -16))),
          'HG8':(((-143,-117),(  73,  97)), (( -69, -51),( -29, -15)),
                 ((-104, -76),(  20,  40)), ((-113, -87),( 108, 136)),
                 ((-104, -76),( 107, 133))),
          'HG9':(((-143,-117),(  73,  97)), (( -69, -51),( 135, 165)),
                 ((-104, -76),( -32, -12)), ((-113, -87),( -40, -14)),
                 ((-104, -76),( -44, -16))),
          'HGA':(((-143,-117),(  73,  97)), (( -69, -51),( 135, 165)),
                 ((-104, -76),( -32, -12)), ((-113, -87),( -40, -14)),
                 ((-104, -76),( 107, 133))),
          'HGB':(((-143,-117),(  73,  97)), (( -69, -51),( 135, 165)),
                 ((-104, -76),( -32, -12)), ((-113, -87),( 108, 136)),
                 ((-104, -76),( -44, -16))),
          'HGC':(((-143,-117),(  73,  97)), (( -69, -51),( 135, 165)),
                 ((-104, -76),( -32, -12)), ((-113, -87),( 108, 136)),
                 ((-104, -76),( 107, 133))),
          'HGD':(((-143,-117),(  73,  97)), (( -69, -51),( 135, 165)),
                 ((-104, -76),(  20,  40)), ((-113, -87),( -40, -14)),
                 ((-104, -76),( -44, -16))),
          'HGE':(((-143,-117),(  73,  97)), (( -69, -51),( 135, 165)),
                 ((-104, -76),(  20,  40)), ((-113, -87),( -40, -14)),
                 ((-104, -76),( 107, 133))),
          'HGF':(((-143,-117),(  73,  97)), (( -69, -51),( 135, 165)),
                 ((-104, -76),(  20,  40)), ((-113, -87),( 108, 136)),
                 ((-104, -76),( -44, -16))),
          'HGG':(((-143,-117),(  73,  97)), (( -69, -51),( 135, 165)),
                 ((-104, -76),(  20,  40)), ((-113, -87),( 108, 136)),
                 ((-104, -76),( 107, 133))),
          'P2A':(((-105, -30),( 100, 180)),),
          'P2B':(((-105, -30),(-180,-160)),),
          'AL_':(((  35, 100),( -25,  60)),),
          'DPA':(((-150,-110),(  45,  90)), ((-135, -30),( -65,  40))),
          'A2D':(((-135, -30),( -65,  40)), ((-135, -30),( -65,  40)),
                 ((-150,-110),(  45,  90))),
          'H__':(((-74,-46),(-54,-26)),),
          '14_':(((-180, 180),(-180, 180)),(( -90, -30),( -60,   0)),
                 ((-120, -60),( -30,  30)),((-180, 180),(-180, 180))),
          '14p':(((-180, 180),(-180, 180)),((  30,  90),(   0,  60)),
                 ((  60, 120),( -30,  30)),((-180, 180),(-180, 180))),
          '24_':(((-180, 180),(-180, 180)),(( -90, -30),(  90, 150)),
                 ((  50, 110),( -30,  30)),((-180, 180),(-180, 180))),
          '24p':(((-180, 180),(-180, 180)),((  30,  90),(-150, -90)),
                 ((-110, -50),( -30,  30)),((-180, 180),(-180, 180))),
          '34_':(((-180, 180),(-180, 180)),(( -90, -30),( -60,   0)),
                 (( -90, -30),( -60,   0)),((-180, 180),(-180, 180))),
          '34p':(((-180, 180),(-180, 180)),((  30,  90),(   0,  60)),
                 ((  30,  90),(   0,  60)),((-180, 180),(-180, 180)))
          }
          
basinBounds  = {'A':(((-135, -30),( -100, 50)),),
                'B':(((-180,-105),(  90, 180)),
                     ((-180,-105),(-180,-170))),
                'L':(((  35, 100),( -25,  60)),),
                'P':(((-105, -30),( 100, 180)),
                     ((-105, -30),(-180,-160))),
                'G':(((-100, -70),(  50, 100)),),         
                'D':(((-150,-110),(  45,  90)),),
                'T':(((  45,  70),(-180,-120)),
                     ((  45,  70),( 160, 180)))}

basins       = ['A','B','P','L','G','D','T']
          
class Coil:

    def __init__(self,f,fname=''):

        self.rdict   = {};
        self.resnums = ();
        self.inserts = ();
        self.resis   = ();
        self.phis    = ();
        self.psis    = ();
        self.omes    = ();
        self.chis    = ();
        self.oldMeso = ();
        self.oldSS   = ();
        self.newMeso = ();
        self.newSS   = ();
        self.name    = fname;
        self.pdbids  = []

        #Set 

        self.resnums,self.inserts,self.resis,self.phis,self.psis,self.omes,\
        self.chis,self.oldMeso,self.oldSS,self.newMeso,self.newSS = \
        zip(*[[x[0:8].strip(),x[9],x[11:14].strip(),float(x[15:25].strip()),
                                         float(x[26:36].strip()),
                                         float(x[37:47].strip()),
                                        (float(x[48:58].strip()),
                                         float(x[59:69].strip()),
                                         float(x[70:80].strip()),
                                         float(x[81:91].strip())),
                                         x[92],x[94],
                                         x[96:98].strip(),x[99]] for x in f \
                                         if len(x.split()) >= 13]);
        

        for i in xrange(len(self.resnums)):
            self.rdict[(self.resnums[i]+self.inserts[i]).strip()] = i;
            self.pdbids.append((self.resnums[i]+self.inserts[i]).strip())

        self.pdbids = tuple(self.pdbids)
            

    #Find all protein dihedrals within a given phi-psi range
    def dihedrals_in_range(self,phimin,phimax,psimin,psimax):

        return [x for x in xrange(len(self.phis)) if \
                phimin<=self.phis[x]<=phimax and \
                psimin<=self.psis[x]<=psimax]

    #Truth statement of whether dihedrals in a particular residue are within
    #a given phi-psi range
    def resi_dihedrals_in_range(self,phimin,phimax,psimin,psimax,idx):

        return (phimin<=self.phis[idx]<=phimax and \
                psimin<=self.psis[idx]<=psimax)

    #Truth statement of whether a given tor file contains a given amino acid
    def has_aa(self,aaName):

        return [x for x in self.resis if x==aaName] != []

    #Returns a formatted string of the phi-psi angles of a given residue
    def phi_psi_str(self,idx):

        return '%10.4f %10.4f' %(self.phis[idx],self.psis[idx])

    #Returns a formatted string of a line in a given .tor entry
    def _string(self,idx):

        return '%8i%1s %3s %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %1s %1s %2s %2s' %(int(self.resnums[idx]),self.inserts[idx],self.resis[idx],
                self.phis[idx],self.psis[idx],self.omes[idx],self.chis[idx][0],
                self.chis[idx][1],self.chis[idx][2],
                self.chis[idx][3],self.oldMeso[idx],self.oldSS[idx],
                self.newMeso[idx],self.newSS[idx])

    #Returns the indices of resides with user-specified names
    def resis_with_names(self,*names):

        return [x for x in xrange(len(self.resis)) if self.resis[x] in names]
    
    #Removes the secondary structure flanking a given coil entry
    def remove_ss(self):

        start_res = 0;
        end_res   = -1;
        
        if self.newSS[start_res] in ['E','H']:
            start_res = 2;

        if self.newSS[end_res] in ['E','H']:    
            end_res  = -2;

        #If there is no secondary structure at the end,
        #include all ending residues
        if end_res == -1:
            loop_list = [(self.resnums[x]+self.inserts[x]).strip() for x in xrange(start_res,len(self.resis))];

        else:
            loop_list = [(self.resnums[x]+self.inserts[x]).strip() for x in xrange(start_res,len(self.resis)-2)]

        return loop_list
    
    #Determines the motif a given residue is in.  Has not been tested extensively.
    def motif(self, rnum, motf=motif_list, pos='f'):

        if self.oldSS[self.rdict[rnum]] == 'H':
            return 'HLX'
        elif self.oldSS[self.rdict[rnum]] == 'E':
            return 'EXT'
        elif self.oldMeso[self.rdict[rnum]] == '*':
            return 'CIS'

        M = motifs

        for i in xrange(len(motf)):

            if (pos == 'f' and \
                self.rdict[rnum]+len(M[motf[i]]) <= len(self.resis)) and \
                in_motif(self.rdict[rnum],self.rdict[rnum]+len(M[motf[i]]),
                         motf[i], self):
                
                return motf[i]

            elif (pos == 'l' and self.rdict[rnum]-len(M[motf[i]]) >= 0) and \
                  in_motif(self.rdict[rnum]-len(M[motf[i]])+1,
                           self.rdict[rnum]+1,motf[i], self):
                return motf[i]

        if (pos == 'l' and self.resis[self.rdict[rnum]] == 'PRO' and \
            self.rdict[rnum]-len(M['A2D']) >= 0 and \
            in_motif(self.rdict[rnum]-len(M['A2D']),
                     self.rdict[rnum],'A2D', self)):
            return 'A2D'


        elif (pos == 'f' and self.resis[self.rdict[rnum]] == 'PRO' and \
              self.rdict[rnum]-len(M['DPA'])+1<= len(self.resis)\
              and in_motif(self.rdict[rnum]-len(M['DPA'])+1,
                           self.rdict[rnum]+1,'DPA', self)):
            return 'DPA'
            
        return 'CL_'

    #Determines the non-secondary structural motif a given residue resides in.
    def noSSMotif(self, rnum, motf=motif_list, pos='f'):

        M = motifs

        for i in xrange(len(motf)):

            if (pos == 'f' and \
                self.rdict[rnum]+len(M[motf[i]]) <= len(self.resis)) and \
                in_noSSMotif(self.rdict[rnum],self.rdict[rnum]+len(M[motf[i]]),
                         motf[i], self):
                
                return motf[i]

            elif (pos == 'l' and self.rdict[rnum]-len(M[motf[i]]) >= 0) and \
                  in_noSSMotif(self.rdict[rnum]-len(M[motf[i]])+1,
                           self.rdict[rnum]+1,motf[i], self):
                return motf[i]

        if (pos == 'l' and self.resis[self.rdict[rnum]] == 'PRO' and \
            self.rdict[rnum]-len(M['A2D']) >= 0 and \
            in_noSSMotif(self.rdict[rnum]-len(M['A2D']),
                     self.rdict[rnum],'A2D', self)):
            return 'A2D'


        elif (pos == 'f' and self.resis[self.rdict[rnum]] == 'PRO' and \
              self.rdict[rnum]-len(M['DPA'])+1<= len(self.resis)\
              and in_noSSMotif(self.rdict[rnum]-len(M['DPA'])+1,
                           self.rdict[rnum]+1,'DPA', self)):
            return 'DPA'
        
        if self.oldSS[self.rdict[rnum]] == 'H':
            return 'HLX'
        elif self.oldSS[self.rdict[rnum]] == 'E':
            return 'EXT'
        elif self.oldMeso[self.rdict[rnum]] == '*':
            return 'CIS'
            
        return 'CL_'

    #Determines the whether the bfactor of a given residue is below a specific
    #threshold.  bb_only=1 means only look at backbone atoms.
    def low_bfac(self,thresh,rnum,bb_only = 0):

        m = pdbparse.read_pdb(coil_path+tor_to_pdb(self.name));

        ch = get_chain(m,get_pdb_chain(self.name))
        ridx = get_resi_idx(ch,rnum)       

        if not bb_only:
            atms = ch[ridx]
        else:
            atms = ch[ridx].atoms_with_name('N','CA','C','O')

        for atm in atms:
            if atm.bf > thresh:
                return 0

        return 1

#Determines whether a given coil segment forms a given motif using its sequence
#of phi-psi angles.
def in_motif(sidx, eidx, mID, coil):

    for i in xrange(sidx,eidx):

        if not in_interval(motifs[mID][i-sidx][0][0],
                           motifs[mID][i-sidx][0][1],
                           coil.phis[i]) or \
           not in_interval(motifs[mID][i-sidx][1][0],
                           motifs[mID][i-sidx][1][1],
                           coil.psis[i]) or \
           coil.oldSS[i] in ['H','E'] or \
           coil.oldMeso[i] == '*':
            return 0

    return 1

#Determines the name of the basin in which the phi-psi angles of a given residue
#reside.
def basinName(phi,psi,rname):

    for bName in basins:
        for bounds in basinBounds[bName]:
            if bounds[0][0]<=phi<=bounds[0][1] and \
               bounds[1][0]<=psi<=bounds[1][1]:
                return bName

    if rname == 'GLY' and phi >= 0:
        return 'Y'

    return 'O'

#Truth statement regarding whether a given residue residues in a given basin.
def inBasin(phi,psi,rname,bName):

    if bName not in ['Y','O']:
        for bounds in basinBounds[bName]:
            if bounds[0][0]<=phi<=bounds[0][1] and \
                   bounds[1][0]<=psi<=bounds[1][1]:
                return 1

        return 0

    else:
        if basinName(phi,psi,rname) == bName:
            return 1

        return 0


#Truth statement regarding whether a given range of residues forms a given motif.
#This function has not been tested extensively.
def in_noSSMotif(sidx,eidx,mID,coil):
    
    for i in xrange(sidx,eidx):

        #print mID, coil.phis[i], coil.psis[i]

        if not in_interval(motifs[mID][i-sidx][0][0],
                           motifs[mID][i-sidx][0][1],
                           coil.phis[i]) or \
           not in_interval(motifs[mID][i-sidx][1][0],
                           motifs[mID][i-sidx][1][1],
                           coil.psis[i]):
            return 0

    return 1
    
def in_interval(min,max,num):

    return min<=num<=max

def coil_length(name):
        return int(string.split(name,'.')[1]);

def get_pdb_id(name):
        return string.split(name,'/')[3][0:4];

def get_pdb_chain(name):
        return string.split(name,'/')[3][4];

def get_resnum(name):    
        return string.strip(string.split(name,'.')[2],'_');


def get_tor_chain(tor_file,chain):    
    if chain == '0' or chain == '_':
        return tor_file
    
    start,end = get_chain_indices(tor_file,chain);
    pdb_chain = -9999;
    if start != -1 and end != -1:
        pdb_chain = tor_file[start:end];

    if pdb_chain == -9999 and chain in ['A','1']:
            return tor_file

    return pdb_chain

def get_chain_indices(tor_file,chain):    
    start = end = -1;
    start_not_found = 1;
    i = 0;
    
    #Get lines for starting and ending the file
    while start_not_found and i < len(tor_file):
        line = tor_file[i].split();
        if len(line) == 2 and line[0] == chain:
            start,end = i+1,int(line[1])+i+1;
            start_not_found = 0;
        i = i+1;
        
    return start,end

#Converts tor id to pdb id
def tor_to_pdb(name):

    name_pieces = string.split(name,'.');

    return name_pieces[0]+'.'+name_pieces[1]+'.'+name_pieces[2]+'.pdb.'+name_pieces[4];

#Gets the array index of a given residue name
def get_resi_idx(chain,rnum):

    for i in xrange(len(chain)):
        if (chain[i].idx+chain[i].icode).strip() == rnum:
            return i

def get_chain(mol,chain):

    if len(mol) == 1:
        return mol[0]

    for i in mol:
        if i.name == chain:
            return i

        

                                                                         

        
        
        
            

        

        
