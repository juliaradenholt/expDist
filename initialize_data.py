import numpy as np
from Bio import AlignIO
from Bio import Align
import Bio.Phylo as ph
import models as m

aa = [0, 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19]

aminoacids = { 'A' : 0, 'R' : 1, 'N' : 2, 'D' : 3, 'C' : 4,
'E' : 5, 'Q' : 6, 'G' : 7, 'H' : 8, 'I' : 9,
'L' : 10, 'K' : 11, 'M' : 12, 'F' : 13, 'P' : 14, 'S' : 15, 'T' : 16,
'W' : 17, 'Y' : 18, 'V' : 19}

#0  A Alanine
#1  R Arginine
#2  N Asparagine
#3  D Aspartic
#4  C Cysteine
#5  E Glutamic
#6  Q Glutamine
#7  G Glycine
#8  H Histidine
#9  I Isoleucine
#10 L Leucine
#11 K Lysine
#12 M Methionine
#13 F Phenylalanine
#14 P Proline
#15 S Serine
#16 T Threonine
#17 W Tryptophan
#18 Y Tyrosine
#19 V Valine

#sample of distances

aligner = Align.PairwiseAligner()

WAG = m.WAG()
LG = m.LG()
PAM =  m.Dayhoff()
BLOSUM = m.BLOSUM62()
JTT = m.JTT()
CPREV = m.cpREV()
MTART = m.mtArt()
MTREV = m.MtREV()

sub_models = [WAG, LG, PAM, JTT]

def readFasta(filename, model):
    """ Reads a file with FASTA format and returns
        the position-differ matrix and the poisson correction distance"""
    align = AlignIO.read(filename, "fasta")
    max_c = [(p1,p2) for p1 in range(len(align)) for p2 in range(p1+1,len(align))]
    names = [(align[i[0]].name, align[i[1]].name) for i in max_c]

    return align, names, max_c

def initializeSequenceData(s1, s2, model):
    """ Takes two sequences and returns a matrix
        corresponds to where position differ and
        the poissoncorrection distance. """

    aa = np.zeros((20,20)) # pairwise alignment scoring
    mutations = 0
    equilibrium_freq = []
    for acid in range(len(s1)):
        i = aminoacids.get(s1[acid])
        j = aminoacids.get(s2[acid])
        equilibrium_freq.append(model.freq[i])
        aa[i][j] += 1
        if i != j:
            mutations += 1
    pdistance = mutations/len(s1)
    pc = poissonCorrection(pdistance)
#    print("Poisson correction distance:", pc)
    return aa, dist(pc), equilibrium_freq

def poissonCorrection(p):
    """ The function takes the p-distance and
        returns the poission correction distance -ln(1-p)"""
    return -1*np.log(1-p)

def dist(pc):
    """ Creates a set of distances from a poisson correction distance value."""

    #distances = [x for x in np.arange(0, 3.0, 0.008)] #uncomment for slow estimator
    if pc < 0.1:
        distances = [x for x in np.arange(0, 0.2, 0.008)]
    elif pc < 0.4:
        distances = [x for x in np.arange(0.9*pc, pc*1.4, 0.015)]
    elif pc > 1.35:
        distances = [x for x in np.arange(pc*0.8, 3, 0.4)]
    else:
        distances = [x for x in np.arange(pc*0.85, pc*1.8, 0.15)]

#    print(len(distances))


    return distances

def draw_tree(filename):
    """ Takes a filename (file in newick format) as input and returns
        a plot of the given tree with matplotlib """
    tree = ph.read(filename, "newick")
    ph.draw(tree,  branch_labels=lambda c: c.branch_length)
