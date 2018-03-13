#! /misc/local/python-2.7.11/bin/python

import sys, getopt, gzip, os
#import pylab as plt
from LINUS import evaluator, molecule
from biomol import pdbparse
from biomol.pdbout import pack_pdb_line
from coil_lib_v3 import *
from math import fabs

#Use the command below if you want to reroute the torsion path from its
#default (/databses/torsions/pdb):
#tor_path = ''


code_3to1 = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q',
             'GLU':'E','GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K',
             'MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W',
             'TYR':'Y','VAL':'V','ASX':'X','GLX':'y'};

'''DGtr of sidechains and the peptide unit from water to 1 molar osmolytes in
   calories per mole per molar.'''
OsmoLib = {
    'A': {'TMAO': -14.64, 'Sarcosine':  10.91, 'Betaine':   4.77,
          'Proline':  -0.07, 'Glycerol':   7.77, 'Sorbitol':  16.57,
          'Sucrose':  22.05,  'Trehalose':  33.25,  'Urea':   0.63},
    'C': {'TMAO':   0.,   'Sarcosine':   0.,   'Betaine':   0.,
          'Proline':   0.,   'Glycerol':   0.,   'Sorbitol':   0.,
          'Sucrose':   0.,    'Trehalose':   0.,    'Urea':   0.},
    'D': {'TMAO': -66.67, 'Sarcosine': -14.20, 'Betaine':-116.56,
          'Proline': -90.51, 'Glycerol': -85.45, 'Sorbitol': -83.88,
          'Sucrose': -37.17,  'Trehalose': -96.54,  'Urea':  43.82},
    'E': {'TMAO': -83.25, 'Sarcosine': -12.61, 'Betaine':-112.08,
          'Proline': -89.17, 'Glycerol': -74.20, 'Sorbitol': -70.05,
          'Sucrose': -41.65,  'Trehalose': -85.92,  'Urea':  40.89},
    'F': {'TMAO':  -9.32, 'Sarcosine': -12.64, 'Betaine':-112.93,
          'Proline': -71.26, 'Glycerol':  59.78, 'Sorbitol':  26.38,
          'Sucrose': -96.35,  'Trehalose': -17.88,  'Urea': -42.84},
    'G': {'TMAO':   0.,   'Sarcosine':   0.,   'Betaine':   0.,
          'Proline':   0.,   'Glycerol':   0.,   'Sorbitol':   0.,
          'Sucrose':   0.,    'Trehalose':   0.,    'Urea':   0.},
    'H': {'TMAO':  42.07, 'Sarcosine': -20.80, 'Betaine': -35.97,
          'Proline': -45.10, 'Glycerol': -17.16, 'Sorbitol': -42.45,
          'Sucrose':-118.66,  'Trehalose': -98.75,  'Urea': -10.24},
    'I': {'TMAO': -25.43, 'Sarcosine':  39.98, 'Betaine':  -1.27,
          'Proline':  -2.72, 'Glycerol':  36.23, 'Sorbitol':  36.90,
          'Sucrose':  28.12,  'Trehalose':  79.66,  'Urea':   1.84},
    'K': {'TMAO':-110.23, 'Sarcosine': -27.42, 'Betaine':-171.99,
          'Proline': -59.87, 'Glycerol': -34.00, 'Sorbitol': -32.47,
          'Sucrose': -39.60,  'Trehalose': -50.08,  'Urea':  17.51},
    'L': {'TMAO':  11.62, 'Sarcosine':  38.33, 'Betaine': -17.73,
          'Proline':   4.77, 'Glycerol': -34.41, 'Sorbitol':  39.07,
          'Sucrose':  37.11,  'Trehalose':  96.17,  'Urea': -14.30},
    'M': {'TMAO':  -7.65, 'Sarcosine':   8.18, 'Betaine': -14.16,
          'Proline': -35.12, 'Glycerol':  13.88, 'Sorbitol':  20.97,
          'Sucrose':  -6.66,  'Trehalose':  29.19,  'Urea':  -8.07},
    'N': {'TMAO':  55.69, 'Sarcosine': -40.93, 'Betaine':  33.17,
          'Proline': -17.71, 'Glycerol':  51.57, 'Sorbitol': -21.21,
          'Sucrose': -28.28,  'Trehalose':  48.67,  'Urea':   1.48},
    'P': {'TMAO':-137.73, 'Sarcosine': -34.23, 'Betaine':-125.16,
          'Proline': -63.96, 'Glycerol': -60.54, 'Sorbitol':  -4.48,
          'Sucrose': -73.02,  'Trehalose': -94.67,  'Urea':  22.62},
    'Q': {'TMAO':  41.41, 'Sarcosine': -10.19, 'Betaine':   7.57,
          'Proline': -32.26, 'Glycerol':  -2.75, 'Sorbitol': -23.98,
          'Sucrose': -40.87,  'Trehalose': -36.34,  'Urea': -14.54},
    'R': {'TMAO':-109.27, 'Sarcosine': -32.24, 'Betaine':-109.45,
          'Proline': -60.18, 'Glycerol': -30.74, 'Sorbitol': -24.65,
          'Sucrose': -79.32,  'Trehalose': -50.33,  'Urea':  19.10},
    'S': {'TMAO': -39.04, 'Sarcosine': -27.98, 'Betaine': -41.85,
          'Proline': -33.49, 'Glycerol':   6.31, 'Sorbitol':  -1.58,
          'Sucrose':  -2.79,  'Trehalose':  -0.98,  'Urea':  19.71},
    'T': {'TMAO':   3.57, 'Sarcosine':  -7.54, 'Betaine':  -0.33,
          'Proline': -18.33, 'Glycerol':  17.53, 'Sorbitol':  13.20,
          'Sucrose':  20.82,  'Trehalose':  26.32,  'Urea':  18.18},
    'V': {'TMAO':  -1.02, 'Sarcosine':  29.32, 'Betaine': -19.63,
          'Proline':   7.96, 'Glycerol':  -1.37, 'Sorbitol':  24.65,
          'Sucrose':  33.92,  'Trehalose':  96.79,  'Urea':  18.62},
    'W': {'TMAO':-152.88, 'Sarcosine':-113.03, 'Betaine':-369.93,
          'Proline':-198.37, 'Glycerol':-122.65, 'Sorbitol': -67.23,
          'Sucrose':-215.27,  'Trehalose':-206.30,  'Urea':-101.19},
    'Y': {'TMAO':-114.32, 'Sarcosine': -26.37, 'Betaine':-213.09,
          'Proline':-138.41, 'Glycerol':-149.50, 'Sorbitol': -53.50,
          'Sucrose': -78.41,  'Trehalose': -80.32,  'Urea':  -4.81},
}

'''From DKP model'''
BackboneLib = {'TMAO': 90., 'Sarcosine':  52., 'Betaine':   67.,
               'Proline':  48., 'Glycerol':   14., 'Sorbitol':  35.,
               'Sucrose':  62.,  'Trehalose':  62.,  'Urea':   -39.}

'''from Creamer et al, Biochemistry, 36:2832-2835 (1997), Table 1'''
unfolded_asaLib = {
    'A': {'BBUpper': 35.9,'BBLower': 19.8, 'SCUpper':  63.6, 'SCLower':  46.6},
    'C': {'BBUpper': 34.5,'BBLower': 18.2, 'SCUpper':  83.0, 'SCLower':  62.9},
    'D': {'BBUpper': 33.9,'BBLower': 18.1, 'SCUpper':  94.8, 'SCLower':  79.2},
    'E': {'BBUpper': 33.5,'BBLower': 17.9, 'SCUpper': 123.9, 'SCLower': 102.8},
    'F': {'BBUpper': 33.3,'BBLower': 15.3, 'SCUpper': 139.8, 'SCLower': 118.7},
    'G': {'BBUpper': 75.7,'BBLower': 54.6, 'SCUpper':   0.0, 'SCLower':   0.0},
    'H': {'BBUpper': 33.4,'BBLower': 14.9, 'SCUpper': 119.1, 'SCLower': 103.9},
    'I': {'BBUpper': 24.7,'BBLower': 15.2, 'SCUpper': 134.1, 'SCLower': 100.1},
    'K': {'BBUpper': 33.8,'BBLower': 18.3, 'SCUpper': 158.8, 'SCLower': 142.5},
    'L': {'BBUpper': 30.7,'BBLower': 14.7, 'SCUpper': 117.7, 'SCLower': 101.4},
    'M': {'BBUpper': 33.8,'BBLower': 16.7, 'SCUpper': 139.5, 'SCLower': 105.3},
    'N': {'BBUpper': 32.7,'BBLower': 17.6, 'SCUpper':  95.6, 'SCLower': 84.5},
    'P': {'BBUpper': 26.1,'BBLower': 18.9, 'SCUpper':  90.5, 'SCLower':  83.5},
    'Q': {'BBUpper': 33.4,'BBLower': 17.2, 'SCUpper': 128.7, 'SCLower': 105.0},
    'R': {'BBUpper': 33.0,'BBLower': 17.1, 'SCUpper': 185.3, 'SCLower': 156.9},
    'S': {'BBUpper': 35.0,'BBLower': 23.8, 'SCUpper':  73.3, 'SCLower':  59.7},
    'T': {'BBUpper': 29.5,'BBLower': 18.6, 'SCUpper':  91.2, 'SCLower':  77.3},
    'V': {'BBUpper': 24.9,'BBLower': 15.9, 'SCUpper': 110.9, 'SCLower':  81.8},
    'W': {'BBUpper': 32.0,'BBLower': 15.1, 'SCUpper': 158.4, 'SCLower': 154.7},
    'Y': {'BBUpper': 33.5,'BBLower': 17.7, 'SCUpper': 152.3, 'SCLower': 131.0}
}    

"""From Creamer and Rose."""
glyXgly_asaLib = {
    'A': {'BB': 46.2, 'SC':  71.9},
    'C': {'BB': 42.6, 'SC': 103.5},
    'Q': {'BB': 37.8, 'SC': 155.4},
    'D': {'BB': 40.5, 'SC': 118.2},
    'E': {'BB': 37.8, 'SC': 148.4},
    'F': {'BB': 38.4, 'SC': 184.4},
    'G': {'BB': 88.1, 'SC':   0.0},
    'H': {'BB': 40.4, 'SC': 162.1},
    'I': {'BB': 30.9, 'SC': 150.1},
    'K': {'BB': 38.7, 'SC': 187.1},
    'L': {'BB': 35.3, 'SC': 157.8},
    'M': {'BB': 38.6, 'SC': 164.8},
    'N': {'BB': 40.2, 'SC': 125.3},
    'P': {'BB': 35.6, 'SC': 111.0},
    'Q': {'BB': 37.8, 'SC': 155.4},
    'R': {'BB': 39.1, 'SC': 216.9},
    'S': {'BB': 44.0, 'SC':  85.8},
    'T': {'BB': 37.9, 'SC': 114.6},
    'V': {'BB': 36.1, 'SC': 128.4},
    'W': {'BB': 37.4, 'SC': 228.9},
    'Y': {'BB': 38.7, 'SC': 198.1}}

atfmt = 'ATOM  %5i %4s%1s%3s %1s%4s%1s   %8.3f%8.3f%8.3f%6.2f%6.2f\n'

class MValueEvaluator:

    def __init__(self,cosolvents):

        self.folded_BB_ASAs = []
        self.folded_SC_ASAs = []

        self.delta_BB_ASAs = []
        self.delta_SC_ASAs = []

        self.mvalues = {}
        self.bbmvals = {}
        self.scmvals = {}

        self.cosolvents = cosolvents

        for c in self.cosolvents:

            self.mvalues[c] = []
            self.bbmvals[c] = []
            self.scmvals[c] = []

        self.ss_assignments = None
        self.continuous_ss_segments = None

    def residue_ASA(self,m,startidx,endidx,asa_ev,rescode):

        bbAtmNams = [' N  ',' CA ',' C  ',' O  ',' H  ']

        folded_BB_ASA = 0.0
        folded_SC_ASA = 0.0
        
        for j in xrange(startidx,endidx):

            if m.atmnam[j] in bbAtmNams:
                folded_BB_ASA += asa_ev.atom(j)
                
            else:
                folded_SC_ASA += asa_ev.atom(j)

        self.folded_BB_ASAs.append(folded_BB_ASA)
        self.folded_SC_ASAs.append(folded_SC_ASA)

    def residue_dgXfer(self,rescode,cosolvent,idx):

        sc_gxe = 0.0
        bb_gxe = 0.0

        unfolded_BB_ASA = average(unfolded_asaLib[rescode]['BBUpper'],
                                  unfolded_asaLib[rescode]['BBLower'])
        unfolded_SC_ASA = average(unfolded_asaLib[rescode]['SCUpper'],
                                  unfolded_asaLib[rescode]['SCLower'])

        self.delta_BB_ASAs.append(unfolded_BB_ASA-self.folded_BB_ASAs[idx])
        self.delta_SC_ASAs.append(unfolded_SC_ASA-self.folded_SC_ASAs[idx])

            
        bb_gxe = BackboneLib[cosolvent]*(unfolded_BB_ASA-\
                                          self.folded_BB_ASAs[idx])/\
                                          glyXgly_asaLib[rescode]['BB']

        if rescode != 'G':        
            sc_gxe = OsmoLib[rescode][cosolvent]*\
                     (unfolded_SC_ASA-self.folded_SC_ASAs[idx])/\
                     glyXgly_asaLib[rescode]['SC']

        self.mvalues[cosolvent].append(bb_gxe+sc_gxe)
        self.bbmvals[cosolvent].append(bb_gxe)
        self.scmvals[cosolvent].append(sc_gxe)

    def smooth_mValues(self,hwindow,cosolvent):

        mVs = self.mvalues[cosolvent]

        return [sum(mVs[x-(hwindow):x+(hwindow)+1])/\
                (hwindow*2+1) for x in xrange(hwindow,len(mVs)-hwindow)]

    def pdbWithMValues(self,m,hwindow,cosolvent,id,mvals):

        resis = [m.resnum[x] for x in m.resptr[:-1]]
        mvalues_avg = [0.0]*hwindow+mvals+\
                      [0.0]*hwindow
        bfact = []

        for i in xrange(len(m.resptr)-1):
            bfact += [mvalues_avg[i]]*(m.resptr[i+1]-m.resptr[i])

        str_out = '%s_%s_%i.pdb' %(id,cosolvent,2*hwindow+1)

        m.writePDB(str_out,bf=bfact)

    def ss_assign(self,pdb_name,resis2assign,hwindow):
        #try:  
            #cf = gzip.open(tor_path+pdb_name[:-1]+'.tor.gz').read().splitlines()
        #except:
            #self.ss_assignments = None
            #return
        
        c = Coil(get_tor_chain(cf,pdb_name[-1]))

        self.ss_assignments = {}
        self.ss_assignments['helix'] = []
        self.ss_assignments['strand']= []
        self.ss_assignments['coil']  = []

        for i in xrange(hwindow,len(resis2assign)-hwindow):
        
            if c.newSS[c.rdict[resis2assign[i]]] == 'H':
                self.ss_assignments['helix'].append(i)
            elif c.newSS[c.rdict[resis2assign[i]]] == 'E':
                self.ss_assignments['strand'].append(i)
            else:
                self.ss_assignments['coil'].append(i)

    def _continuous_segments(self,residues,minres):

        self.continuous_ss_segments = []

        segments = self.continuous_ss_segments

        for i in xrange(len(residues)):

            if residues[i] == minres:
                segments.append([residues[i]]);

            elif residues[i]-residues[i-1] != 1:
                segments.append([residues[i]-1,residues[i]])

            else:
                segments[-1].append(residues[i])

    def ssPlot(self,resis,hwindow,mvalu_avgs,symbol):

        leg = None

        if not self.ss_assignments:
            plt.hold(True)
            leg, = plt.plot(resis[hwindow:-hwindow],mvalu_avgs)
            return leg

        for i in ['helix','strand','coil']:

            self._continuous_segments(self.ss_assignments[i],hwindow)

            for j in xrange(len(self.continuous_ss_segments)):
            
                curResis = [resis[x] for x in self.continuous_ss_segments[j]]
                try:
                    curmvals = [mvalu_avgs[x-hwindow] for x in \
                                self.continuous_ss_segments[j]]
                except IndexError:
                    offset = self.continuous_ss_segments[j][-1]-2-\
                             len(mvalu_avgs)
                    curmvals = [mvalu_avgs[x-hwindow] for x in \
                                self.continuous_ss_segments[j][:-offset]]
                    curResis = curResis[:-offset]

                plt.hold(True)

                if i == 'helix':
                    plt.plot(curResis,curmvals,'m'+symbol)
                    
                elif i == 'strand':
                    plt.plot(curResis,curmvals,'c'+symbol)

                elif i == 'coil':
                    if j == 0:
                        leg, = plt.plot(curResis,curmvals,'k'+symbol)
                    else:
                        plt.plot(curResis,curmvals,'k'+symbol)

        return leg

def usage():

    sys.exit("fragment_plots.py [-h] [options] <pdb or list of pdbs>")

def help():

    sys.exit(\
        """
           Name: fragment_plots.py
           Author: Lauren Porter (lperski2@jhu.edu)
           Date: 1/10/11
           Required packages: Matplotlib, LINUS, biomol, coil_lib_v3, Numeric
           Usage: fragment_plots.py [options] <input file or pdb_ID+chain>

           Description:  This program calculates m-values for whole protein
           chains and the central residues of n-residue segments and plots
           them on a graph.  Additional options are discussed below

           
           options:

           -C: cosolvents.  Full names must be entered as follows:
                  o Urea
                  o TMAO
                  o Trehalose
                  o Sucrose
                  o Proline
                  o Betaine
                  o Sorbitol
                  o Glycerol
                  o Sarcosine

               Default cosolvent is Urea only.
                  
               ***Note: the -C option must be entered for each cosolvent.  IE,
               if m-values are desired for both Urea and Trehalose, one must
               enter -C Urea -C Trehalose.  See examples below.
                  Also: graphs are made for each cosolvent separately to avoid
               cluttered figures.***
               
                

           -s: starting residue.  Use only if you want to get the m-value of a
               protein segment (as opposed to an entire protein chain).
               Argument should be a residue number in the PDB.  Default is None.
               ***Note, this routine has not been well-developed***

           -e: ending residue.  Use only if you want to get the m-value of a
               protein segment (as opposed to an entire protein chain).
               Argument should be a residue number in the PDB.  Default is None
               ***Note, this routine has not been well-developed***

           -1: get the mvalue of just 1 pdb.  0 for no; 1 for yes.  Last
               argument should be the pdb id + chain  (ie 3bdcA for chain A of
               pdb 3bdc).  Note _ is used for PDBs with out chain names.
               Default is 0.  This means inputs should be a text
               file with 1 pdb name+chain on each line.  For example, the file
               could contain:

                   8lyz_
                   3bdcA
                   1stn_

               The program would then calculate m-values and output whatever
               the user specifies for each pdb.
                   

           -U: Use an already clean pdb in current directory.  1 for yes, 0
               for no.

           -W: Window size used for averaging (must be an integer).
               Default is 7.   

           -x: Extension for pdb files.  Those in the Roselab database are
               '.ent.gz' (to save space).  You may need to change the default
               extension depending on your pdb database to '.pdb', for example.

           -p: Path to the pdb file database.  The Roselab path is
               '/databases/pdb/' and an additional 'pdb' is added since all
               pdb IDs are prefixed with it.  You may need to change the
               default path depending on your pdb database.

           -O: Write PDBs for each cosolvent with residue-specific m-values
               in the b-factor column.  1 for yes; 0 for no.  Default is 0.

           -m: Print m-values for each cosolvent.  1 for yes; 0 for no.
               Default is 1.

           -g: Write a figure of dG-xfer versus residue number for all
               cosolvents specified.  Input is the desired name of figure.
               No extension necessary: .png is automatically added.  Secondary
               structure is classified and colored (Magenta, helix;
               cyan, strand; black, coil) if the PROSST secondary structure 
               classifications are installed in /databases/torsions/
               Default:'', meaning no figure will be written.

           -r: Write a text file containing residue-specific deltaG transfers 
               for ***fragments***.  Residue number in the first column, dG in
               the second column, and delta ASAs (total, BB and SC) in the 
               next 3 respective columns. 1 for yes; 0 for no.  Default: 0.
               File is named: {pdb_name}_{cosolvent}_fragMvals_win{window}.txt

           Examples:

               fragment_plots.py -C Urea -C TMAO -g 8lyz_urea_tmao -1 1 -W 9
               8lyz_

               The command above calculates and prints the total m-value of
               the PDB 8lyz (lysozyme) in both urea and TMAO and makes a figure
               of dG xfer for both urea and TMAO versus residue number
               (averaged over a window of 9 residues, residue to be plotted is
               central).

               fragment_plots.py -O 8lyz_urea -1 1 -U 1 8lyz_

               The command above calculates and prints the total m-value of
               the PDB 8lyz (lysozyme) in urea and writes a PDB file with
               residue-specific delta G transfers in the b-factor column.

               fragment_plots.py -C Urea -C TMAO -C Trehalose -m 0
                                  -g test_ pdbs

               The command above calculates the total-mvalues of each protein
               in the file named 'pdbs'.  It then outputs figures of dG xfer
               for urea, TMAO, and trehalose versus residue number (averaged
               over a window of 7 residues, residue to be plotted is central) 
               
                  """)

def get_opts(argv):

    default_opts = {'-C':['Urea',],
                    '-f':'',
                    '-s':None,
                    '-e':None,
                    '-1':0,
                    '-U':0,
                    '-W':7,
                    '-x':'.ent.gz',
                    '-p':'/databases/pdb/pdb',
                    '-O':'',
                    '-m':1,
                    '-g':'',
                    '-r':0}

    default_dos = {'-C':str,
                   '-f':str,
                   '-s':str,
                   '-e':str,
                   '-1':int,
                   '-U':int,
                   '-W':int,
                   '-x':str,
                   '-p':str,
                   '-O':str,
                   '-m':int,
                   '-g':str,
                   '-r':int}

    short_opts = 'C:f:s:e:1:U:P:W:x:p:T:O:m:g:huA:r:'

    if not argv:
        usage()

    opts,args = getopt.getopt(argv,short_opts)

    cosolvents = []

    for o,v in opts:
        if o == '-C':
            cosolvents.append(v)
        elif o=='-h':
            help()
           
        default_opts[o] = default_dos[o](v)

    default_opts['-f'] = args[0]

    if cosolvents:
        default_opts['-C'] = cosolvents

    return default_opts

def get_protein_chain(pdbName,chnName,pdbPath,extension):
    
    try:
        mol = pdbparse.read_pdb(pdbPath+pdbName+extension);

    except IOError:
        return []

    if chnName == '_':
        chain = mol[0]
    else:
        chain = get_chain(mol,chnName);

    return chain

def get_chain(mol,chain):

    for i in mol:
        if i.name == chain:
            return i

def clean_chain(chn):

    chn.delete_hetero()
    chn.delete_alt_locs()
    chn.delete_atoms_with_name('OXT')

def writePDB(mOut, fName):

    fo = open(fName+'.pdb','w')

    for j in xrange(len(mOut)):
        for k in xrange(len(mOut[j])):
            if 'H' not in mOut[j][k].name or mOut[j][k].name in \
                   ['OH','CH3','CH2']:
                pack_pdb_line(mOut[j].atoms()[k],mOut[j][k].idx,
                              mOut[j].name,mOut[j].idx,mOut[j].icode,
                              fName[-1],fo,atfmt)

def get_idx_converter(mol,segment):

    pdbi2listi = {}

    for i in xrange(mol.numres):

        pdbi2listi[(str(mol.pdbidx[i])+mol.pdbicd[i]).strip()] = i

    return pdbi2listi

def average(a,b):

    return (a+b)/2.0

def get_resi_breaks(pdb):

    resi_dict = {}
    oldres = None
    ridx = 0

    for i in xrange(len(pdb)):
        if pdb[i][22:27] != oldres:
            oldres = pdb[i][22:27]
            resi_dict[ridx] = i
            ridx += 1

    resi_dict[ridx] = i+1

    return resi_dict

        
def get_total_mvalues(pdb,cosolvents,offsets):

    bbAtmNams = [' N  ',' CA ',' C  ',' O  ',' H  ']

    #print pdb
    
    saveout = sys.stdout
    sys.stdout = open(os.devnull, 'w')
    m = molecule.MoleculeWithBackboneHydrogens(pdb)
    #m.writePDB('test.pdb')
    sys.stdout = saveout

    #sys.exit()

    asa_ev = evaluator.SurfaceAreaEvaluator(m, m.numres, probe=1.4,
                                            ndiv=3, hphil=1)
    asa = asa_ev.evaluate()

    segment = [offsets[0],m.numres-offsets[1]]
    #print segment

    totalMvalues = MValueEvaluator(cosolvents)
    
    for i in xrange(segment[0],segment[1]):

        #print m.resnam[i]

        startidx = m.resptr[i]

        if i != m.numres:
            endidx = m.resptr[i+1]
        else:
            endidx = m.numatm
            
        try:
            rescode = code_3to1[m.resnam[i]]
        except:
            print 'Warning: skipping unconventional residue %s' %(m.resnam[i])
            continue

        totalMvalues.residue_ASA(m,startidx,endidx,asa_ev,rescode)

        for c in cosolvents:
            totalMvalues.residue_dgXfer(rescode,c,-1)
            
    return totalMvalues

def run(argv):

    opts = get_opts(argv)

    cosolvents = opts['-C']
    files = opts['-f']
    single_pdb = opts['-1']
    use_clean_pdb = opts['-U']
    segStartIdx = opts['-s']
    segEndIdx   = opts['-e']
    pdbPath = opts['-p']
    pdbExtension = opts['-x']
    window = opts['-W']
    writePdbs = opts['-O']
    printTotalMvalues = opts['-m']
    figName = opts['-g']
    residue_mVals = opts['-r']
 
    
    if single_pdb:
        fnames = [files]
    else:
        fnames = open(files).read().splitlines()

    for i in xrange(len(fnames)):

        #print fnames[i]
        #sys.exit()

        if not use_clean_pdb:
        
            chn = get_protein_chain(fnames[i][0:-1],fnames[i][-1],
                                    pdbPath,pdbExtension)
            #print 'Here'
            #print chn
            try:
                clean_chain(chn)
            except:
                'Protein chain %s not found.  Skipping...' %(fnames[i])
                continue

            writePDB(chn,fnames[i])
            del chn
        
        if segStartIdx and segEndIdx:
            segment = [segStartIdx,segEndIdx]
        else:
            segment = []

        if not use_clean_pdb:
            pdbin = open(fnames[i]+'.pdb').read().splitlines()
        else:
            pdbin = open(fnames[i]+'.pdb').read().splitlines()
        
        if segment:
            saveout = sys.stdout
            sys.stdout = open(os.devnull, 'w')
            m = molecule.MoleculeWithBackboneHydrogens(fnames[i]+'.pdb')
            sys.stdout = saveout
            pdbi2listi = get_idx_converter(m,segment)
            segment = [pdbi2listi[x] for x in segment]
            segment[1]+=1
            resi_dict = get_resi_breaks(pdbin)
            offset = 1
            ct_offset = 1
            nt_offset = 1

            if segment[0] == 0:
                offset = 0
            if segment[1] == len(pdbin):
                ct_offset = 0

            if pdbin[segment[0]-offset].split()[3] != 'ACE':
                
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
                #ct_offset=1
                #print segment[1]
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
            
                segpdb =ACE+\
                    pdbin[resi_dict[segment[0]]:resi_dict[segment[1]]]+NME

            segmentMvalues = get_total_mvalues(segpdb,cosolvents,
                                               [1,1])

        else:
                    #print resi_dict.keys()
                    #sys.exit()
            print 'Better.'
            resi_dict = get_resi_breaks(pdbin)
            offset = 0
            ct_offset = -1
            saveout = sys.stdout
            sys.stdout = open(os.devnull, 'w')
            m = molecule.MoleculeWithBackboneHydrogens(fnames[i]+'.pdb')
            sys.stdout = saveout
            segpdb = pdbin[0:]
            
            segmentMvalues = get_total_mvalues(pdbin,cosolvents,[1,1])
            nt_offset = 1
            ct_offset = 1
        
        pdbMvalues = get_total_mvalues(pdbin,cosolvents,\
                                       [pdbi2listi[segStartIdx],\
                                        m.numres-pdbi2listi[segEndIdx]-1])
        print 'To here'


        if printTotalMvalues:
            print 'Yes :)'
            for c in cosolvents:
                segSum = sum(segmentMvalues.mvalues[c])
                situSum= sum(pdbMvalues.mvalues[c])
                print segSum, situSum
                #sys.exit()
                BBDiff = -sum(pdbMvalues.folded_BB_ASAs)+\
                             sum(segmentMvalues.folded_BB_ASAs)
                SCDiff = -sum(pdbMvalues.folded_SC_ASAs)+\
                             sum(segmentMvalues.folded_SC_ASAs)

                #print BBDiff, SCDiff
                mval_ratio = segSum/situSum
                # print 'Resis%6s -%5s %10s m-value: %10.4f %10.4f %10.4f'\
                    #%(opts['-s'],opts['-e'],c,segSum,situSum,mval_ratio)
                #print '%10.4f %10.4f %10.4f' %(BBDiff,SCDiff,BBDiff/(BBDiff+SCDiff))
                bbsegMvals = sum(segmentMvalues.bbmvals[c])
                scsegMvals = sum(segmentMvalues.scmvals[c])
                bbpdbMvals = sum(pdbMvalues.bbmvals[c])
                scpdbMvals = sum(pdbMvalues.scmvals[c])

                for i in xrange(len(segmentMvalues.bbmvals[c])):
                    print segmentMvalues.bbmvals[c][i], pdbMvalues.bbmvals[c][i], segmentMvalues.scmvals[c][i], pdbMvalues.scmvals[c][i]

                #print '%10.4f %10.4f %10.4f' %(bbsegMvals,scsegMvals,
                 #                              fabs(bbsegMvals)/(fabs(bbsegMvals)+fabs(scsegMvals)))
                #print '%10.4f %10.4f %10.4f' %(bbpdbMvals,scpdbMvals,
                 #                              fabs(bbpdbMvals)/(fabs(bbpdbMvals)+fabs(scpdbMvals)))
                #print '%10.4f %10.4f' %(bbsegMvals/bbpdbMvals,
                 #                       scsegMvals/scpdbMvals)

                print sum(pdbMvalues.folded_BB_ASAs),sum(segmentMvalues.folded_BB_ASAs), sum(pdbMvalues.folded_SC_ASAs), sum(segmentMvalues.folded_SC_ASAs)
                print ((bbsegMvals+scsegMvals)/(bbpdbMvals+scpdbMvals))
                print bbsegMvals, bbpdbMvals, scsegMvals, scpdbMvals
                print bbpdbMvals+scpdbMvals

                return bbpdbMvals,scpdbMvals,bbsegMvals,scsegMvals,sum(pdbMvalues.folded_BB_ASAs),sum(pdbMvalues.folded_SC_ASAs), sum(segmentMvalues.folded_BB_ASAs), sum(segmentMvalues.folded_SC_ASAs)
                

        #if not use_clean_pdb:
            #os.remove(fnames[i]+'.pdb')

                    
if __name__ == '__main__':

    run(sys.argv[1:])

    


