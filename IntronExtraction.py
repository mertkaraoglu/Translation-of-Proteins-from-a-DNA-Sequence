# This script extracts the intron from a DNA sequence encoded in the FASTA format
# The link to the file is: https://www.ncbi.nlm.nih.gov/nuccore/NC_000017.11?report=fasta&from=7668402&to=7687550
# The intuition is as follows:
# 1. Find the initial point "ATG"
# 2. Find the Branch Points
# 3. Locate the 5' and 3' intron sites (GT - AG)
# 4. Extract the exons
# 5. Search for the end point "TGA" and cut the last exon in there.

# Required Libraries
import numpy as np
import sys

#all possible translations from exon to polypeptide chain
d = {"TTT":"Phe", "TTC":"Phe","TTA":"Leu", "TTG":"Leu","TCT":"Ser", "TCC":"Ser"
     ,"TCA":"Ser", "TCG":"Ser","CTT":"Leu", "CTC":"Leu","CTA":"Leu", "CTG":"Leu"
    , "CCT": "Pro", "CCC": "Pro", "CCA": "Pro", "CCG": "Pro"
    , "ACT": "Thr", "ACC": "Thr", "ACA": "Thr", "ACG": "Thr"
    , "GCT": "Ala", "GCC": "Ala", "GCA": "Ala", "GCG": "Ala"
    , "GTT": "Val", "GTC": "Val", "GTA": "Val", "GTG": "Val"
    , "ATT": "Ile", "ATC": "Ile", "ATA": "Ile", "ATG": "Met"
    , "TAT": "Tyr", "TAC": "Tyr", "TAA": "Stop", "TAG": "Stop"
    , "CAT": "His", "CAC": "His", "CAA": "Gln", "CAG": "Gln"
    , "AAT": "Asn", "AAC": "Asn", "AAA": "Lys", "AAG": "Lys"
    , "GAT": "Asp", "GAC": "Asp", "GAA": "Glu", "GAG": "Glu"
    , "TGT": "Cys", "TGC": "Cys", "TGA": "Stop", "TGG": "Trp"
    , "CGT": "Arg", "CGC": "Arg", "CGA": "Arg", "CGG": "Arg"
    , "AGT": "Ser", "AGC": "Ser", "AGA": "Arg", "AGG": "Arg"
    , "GGT": "Gly", "GGC": "Gly", "GGA": "Gly", "GGG": "Gly"}

#translating function
def Translate(exons):
    l = []

    for i in range(0,len(exons)-2,3):
        l.append(d[exons[i:i+3]])
        #this is optional if we want the sequence to be
        #printed without any spaces and separation:
        sys.stdout.write(d[exons[i:i+3]])

    #if we print what the function returns we will get separated view
    #of the amino acids:
    sys.stdout.write('\n')
    return l


# The FASTA file
fileContent = open("sequence.fasta", "r")

# Converting the sequence into a char list
seq = ""
for line in fileContent:
    if line.startswith(">"):
        continue
    else:
        seq = seq + line[:-1]

# Getting the initial point of the pre-mRNA
for ind in range(len(seq) - 2):
    if seq[ind:ind + 3] == "ATG":
        init_ind = ind
        break     
seq = seq[init_ind:]

# Finding the locations of the branch points
branch_points = []
for ind in range(len(seq) - 4):
    if seq[ind:ind + 5] == "CTAAC" or seq[ind:ind + 5] == "CTAAT" or seq[ind:ind + 5] == "CTGAC" or seq[ind:ind + 5] == "CTGAT" or seq[ind:ind + 5] == "CCAAC" or seq[ind:ind + 5] == "CCAAT" or seq[ind:ind + 5] == "CCGAC" or seq[ind:ind + 5] == "CCGAT" or seq[ind:ind + 5] == "TTAAC" or seq[ind:ind + 5] == "TTAAT" or seq[ind:ind + 5] == "TTGAC" or seq[ind:ind + 5] == "TTGAT" or seq[ind:ind + 5] == "TCAAC" or seq[ind:ind + 5] == "TCAAT" or seq[ind:ind + 5] == "TCGAC" or seq[ind:ind + 5] == "TCGAT":
        branch_points.append(ind+3)

# Finding intron start and end points
int_init_ind = np.zeros(len(seq), dtype=int)
int_end_ind = np.zeros(len(seq), dtype=int)

int_init_ind = list(int_init_ind)
int_end_ind = list(int_end_ind)


for ind in range(len(branch_points)):
    for i in range(branch_points[ind], branch_points[ind]-60, -1):
        if seq[i:i+2] == "GT":
            int_init_ind[i] = 1
            break
    for j in range(branch_points[ind], branch_points[ind]+200, 1):
        if seq[j:j+2] == "AG":
            int_end_ind[j+1] = 1
            break

# Locating the exon element indices in the sequence
exon_ind = np.ones(len(seq), dtype=int)
exon_ind = list(exon_ind)

for ind in range(len(seq)):
    if int_init_ind[ind] == 1:
        for ind2 in range(ind,len(seq)):
            if int_end_ind[ind2] == 0:
                exon_ind[ind2] = 0
            else:
                break

exon_loc = []
for i, j in enumerate(exon_ind):
    if j == 1:
        exon_loc.append(i)

exons = [seq[i] for i in exon_loc]
exons = ''.join(exons)

# Finding the pre-mRNA end point and cut the last exons there
end_indices = []
for ind in range(len(exons) - 2):
    if exons[ind:ind + 3] == "TGA":
        end_indices.append(ind + 2)
end_ind = end_indices[-1]

exons = exons[:end_ind+1]

#Translate(exons)
print(Translate(exons))
sys.stdout.flush()

