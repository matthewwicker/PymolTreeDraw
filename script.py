#-------------This can be run with the command “run file.py” from within pymol
#-------------You will need to read a PDB file before running this or at that to the program
from pymol import cmd
from pymol.cgo import *
from pymol.vfont import plain

import sys
import time
import argparse
import os
from os.path import basename
from tetrahedron import *
import numpy as np

# SET UP ALL OF THE PATHS WE MAY NEED
seq_name = "RNA-Test-Seq/2HOJA"
sequence_file = seq_name + ".fasta"
tree_file = seq_name + "_ann.outr"
known_pdb = seq_name + ".pdb"

# LOAD IN THE PDB TO PYMOL
cmd.load( known_pdb )

# LOAD IN A SEQUENCE
fasta = open(sequence_file)
s = fasta.readlines()
sequence = s[1]

# LOAD IN THE THREE TREE
tree = load_three_tree(tree_file, sequence)

# LOAD IN THE INITIAL LINES FROM THE PDB FILE

pdbf = open(known_pdb)
s = pdbf.readlines()

begun_read = False
atomic_lines = []
for i in s:
    broken_line = i.split()
    # We do not want to take any hydrogens in to the structure.
    if(broken_line[0] == "TER" and begun_read):
        break
    if(broken_line[0] == "ATOM" or broken_line[0] == "HETATM"):
        if(broken_line[2][0] == 'H'):
            continue
        # For now, let's not take any alternate confirmations
        if(broken_line[2][-1] == 'B' or broken_line[3] == 'B'):
            continue
        atomic_lines.append(i)
    begun_read = True

atomic_lines.append("ATOM    -1  x    x x xxx      xx.xxx  xx.xxx xxx.xxx  x.xx xx.xx           x \n")
intitial_data = atomic_lines[0].split()

nucleotides = []

coords = []
labels = []

if(len(intitial_data) == 11 or len(intitial_data) == 12):
    chain_id = intitial_data[4]
    nt_num = intitial_data[5]

elif(len(important_lines) == 13):
    chain_id = intitial_data[5]
    nt_num = intitial_data[6]

current_chain_id = chain_id
current_nt_nums = nt_num
current_nt_index = 0
first = True

added = 0
for i in atomic_lines:
    if(first ):
        first = False
        continue
    else:
        
        atom_data = i.split()
        
        if(len(atom_data) == 11 or len(atom_data) == 12):
            chain_id = atom_data[4]
            nt_num = atom_data[5]
            x, y, z = 6, 7, 8
        
        if(len(atom_data) == 13):
            chain_id = atom_data[5]
            nt_num = atom_data[6]
            x, y, z = 7, 8, 9
        
        if(current_nt_nums != nt_num):
            temp = nucleotide()
            temp.number = current_nt_index
            temp.coords = coords
            temp.labels = labels
            temp.nt = sequence[current_nt_index]
            nucleotides.append(temp)
            added += 1
            current_nt_index+=1
            current_nt_nums = nt_num
            coords = []
            labels = []
        
        coords.append( [atom_data[x], atom_data[y], atom_data[z]] )
        labels.append(atom_data[2])



"""
   NOW WE HAVE ALL OF THE TREE INFORMATION AND WE HAVE ALL OF THE NUCLEOTIDE COORDINATES 
   SPECIFIED:
   
   nucleotides - LIST OF NUCLEOTIDE OBJECTS FOR EVERY ITEM IN THE FIRST CHAIN
   tree - LIST OF TETRAHEDRON OBJECTS
   
"""



def draw_in_pdb(num1, num2, nucs, constraint, index):
    print num1, num2
    
    cent1 = nucs[num1-1].coords[0]
    cent2 = nucs[num2-1].coords[0]
    
    print cent1, cent2
    
    if(constraint == "NaN"):
        r1,g1,b1 = 1,1,1 # color (white)
        r2,g2,b2 = 1,1,1 # color (white)
    elif("cWW" in constraint):
        r1,g1,b1 = 1,0,1 # color (purple)
        r2,g2,b2 = 1,0,1 # color (purple)
    elif("-BB-" in constraint and "s" in constraint):
        r1,g1,b1 = 0,0,1 # color (blue)
        r2,g2,b2 = 0,1,0 # color (green)
    elif("-BB-" in constraint):
        r1,g1,b1 = 0,0,1 # color (blue)
        r2,g2,b2 = 0,0,1 # color (blue)
    elif("s" in constraint):
        r1,g1,b1 = 0,1,0 # color (green)
        r2,g2,b2 = 0,1,0 # color (green)
    else:
        r1,g1,b1 = 1,0,0 # color (red)
        r2,g2,b2 = 1,0,0 # color (red)
    
    x1,y1,z1 = float(cent1[0]), float(cent1[1]), float(cent1[2]) # start point
    x2,y2,z2 = float(cent2[0]), float(cent2[1]), float(cent2[2]) # start point
    cmd.load_cgo( [ 9.0, x1, y1, z1, x2, y2, z2, radius, r1, g1, b1, r2, g2, b2 ], "cylinder%s"%(index) )


radius = 0.125

index = 0
for i in tree:
    # Draw the x1 - x2 line
    num1 = i.x_1.number
    num2 = i.x_2.number
    draw_in_pdb(num1, num2, nucleotides, i.e_1_2_constraint, index)
    index += 1

    # Draw the x1 - x3 line
    num1 = i.x_1.number
    num2 = i.x_3.number
    draw_in_pdb(num1, num2, nucleotides, i.e_1_3_constraint, index)
    index += 1

    # Draw the x1 - x4 line
    num1 = i.x_1.number
    num2 = i.x_4.number
    draw_in_pdb(num1, num2, nucleotides, i.e_1_4_constraint, index)
    index += 1

    # Draw the x2 - x3 line
    num1 = i.x_2.number
    num2 = i.x_3.number
    draw_in_pdb(num1, num2, nucleotides, i.e_2_3_constraint, index)
    index += 1

    # Draw the x2 - x4 line
    num1 = i.x_2.number
    num2 = i.x_4.number
    draw_in_pdb(num1, num2, nucleotides, i.e_2_4_constraint, index)
    index += 1

    # Draw the x3 - x4 line
    num1 = i.x_3.number
    num2 = i.x_4.number
    draw_in_pdb(num1, num2, nucleotides, i.e_3_4_constraint, index)
    index += 1

# THIS WILL DRAW A CYLINDER BETWEEN ALL BACKBONE EDGES

#for i in range(len(nucleotides)-1):

#    cent1 = nucleotides[i].coords[0]
#    cent2 = nucleotides[i+2].coords[0]
    
#    print cent1, cent2

#    x1,y1,z1 = float(cent1[0]), float(cent1[1]), float(cent1[2]) # start point
#    r1,g1,b1 = 1,1,1 # color (red)
#    x2,y2,z2 = float(cent2[0]), float(cent2[1]), float(cent2[2]) # start point
#    r2,g2,b2 = 1,1,1 # color (red)
#    radius = 0.5
#    cmd.load_cgo( [ 9.0, x1, y1, z1, x2, y2, z2, radius, r1, g1, b1, r2, g2, b2 ], "cylinder%s"%(labels[i]) )

