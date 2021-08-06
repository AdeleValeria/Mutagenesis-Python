# -*- coding: utf-8 -*-
"""
Adele Valeria
Computational Biology and Systems Biology Mini-Project

Programmed using Pyython 3.8.8
"""

# import packages
import os
import random

# change directory
os.chdir("C://mini/6ZXN")

#-----PERFORM MUTATIONS BASED ON POSITIONS OF MUTATION IN NATURAL VARIANTS----
'''
The positions of natural variants for SARS-CoV-2 Spike protein can be
found on: https://www.uniprot.org/uniprot/P0DTC2

file : a FASTA file containing Spike protein sequence
nk : a dictionary containing the positions of mutation and 
     their corresponding replaced residues in natural variants
    
Outcome
-------
The function will write the protein sequence of each natural variant into file
'''
def mutateSpike(file, nk):
    with open (file, mode = "r+") as f:
        for line in f:
            if line.startswith(">"):
                continue
            else:
                for j, k in nk.items():
                    temp = line
                    for i in range (1,len(line)):
                        if i in k:
                            temp = temp[:i-1]+k[i]+temp[i:]
                    f.write('>'+j+'\n'+temp)
    return

#-----------COUNT THE NUMBER OF MUTATIONS FOR EACH NATURAL VARIANT-----------
'''
nv : a dictionary containing the positions of mutations and
     their corresponding replaced residues in natural variants

Outcome
-------
The function will return a dictionary containing the number of mutations 
for each natural variant
'''
def countMutVariant(nv):  
    mut_count = {}                 
    for j, k in nv.items():
        mut_count[j+'_comparison'] = len(k)
    return mut_count

#------------------PERFORM MUTATIONS ON BURIED RESIDUES-----------------------
'''
file: a FASTA file containing Spike protein sequence
S1: positions of buried residues in S1 region (pos 13 - 685)
S2: positions of buried residues in S2 region (pos 686 - 1273)
nv: a dictionary containing the number of mutations for each natural variant
copy : an int indicating the desired number of replicates
AA : amino acid residue to replace the randomly selected buried residues

Outcome
-------
The function will write the protein sequence of each mutant into file. 

Example
-------
B1351 natural variant has 8 mutations. This function will randomly select 8
positions of buried residues and replace them with AA, then write the protein 
sequence into file for n number of replicates.
'''

def mutateSpikeRandom(file, S1, S2, nv, copy, AA):
    buried = [S1,S2]
    random_mut = {}
    for n in range (copy):
        for i, j in nv.items():
            while (True):
                # mutation in S2 region rarely happens, thus lower probability
                mut_buried = random.choices(buried, weights = (80,20), k=j)
                mut_buried = map(random.choices, mut_buried)
                mut_buried = sum(sorted(list(mut_buried)),[])
                if len(mut_buried) == len(set(mut_buried)):
                    break
            random_mut[str(i)+"_"+str(n)+str(mut_buried)] = mut_buried
    with open(file, mode = "r+") as f:
        for line in f:
            if line.startswith(">"):
                continue
            else:
                for j, k in random_mut.items():
                    temp = line
                    for i in range (1,len(line)):
                        if i in k:
                            temp = temp[:i-1]+AA+temp[i:]
                    f.write('>'+j+'\n'+temp)
    return


# fasta file of monomer protein sequence
file = "6ZXN COMP.fasta" 

# list natural variants and their corresponding mutation positions
natural_variant = {"P1": {18:"F",20:"N",26:"S",138:"Y",190:"S",417:"N",484:"K",
                          501:"Y",614:"G",655:"Y",1027:"I"},
                   "B1351": {18:"F",80:"A",215:"G",417:"N",484:"K",501:"Y",614:
                             "G",701:"V"},
                   "19B/501Y": {18:"F",452:"R",501:"Y",653:"V",655:"Y",796:"Y"},
                   "B11318": {95:"I",484:"K",614:"G",796:"H"},
                   "B1525": {52:"R",484:"K",888:"L"}}
                   
'''
Identify the positions of buried residues using PDBe PISA. 
Spike protein has 2 subunits, and mutation occurs more frequently in subunit S1,
thus I separated the buried residues into 2 lists
'''
buried_S1 = [55,79,89,90,91,93,95,103,107,117,131,189,193,223,241,244,260,276,
             277,299,353,400,402,431,438,497,507,511,512,539,552,575,649]
buried_S2 = [694,816,818,822,876,877,906,923,980,1004,1025,1029,1049,1052,1053,
             1054,1060,1062,1063,1066,1067,1080,1081]

# specify the number of replicates for each number of mutations
copy = 3

# specify the amino acid to replace the randomly selected buried residues
AA = 'A'


mutateSpike(file,natural_variant)
mutations_NV = countMutVariant(natural_variant)
mutateSpikeRandom(file, buried_S1, buried_S2, mutations_NV, copy, AA)


