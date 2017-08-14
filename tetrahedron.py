import sqlite3
import operator
import copy
from math import sqrt
import numpy as np
import random
import re
from numpy import linalg as LA

def getDist(vi, vj):
    dist = sqrt((vi[0] - vj[0])**2 + (vi[1] - vj[1])**2 + (vi[2] - vj[2])**2)
    return dist

class nucleotide:
    def __init__(self):
        self.nt = 'N'
        self.coords = []
        self.labels = []
        self.number = -1
    def show(self):
        print(self.number, self.nt, self.labels[0])
    def isIncomplete(self):
        print("Coords:")
        print(self.coords)
        print("Retval:")
        retval = np.isnan(np.asarray(self.coords).astype('float')).any()
        print(retval)
        return retval

class candidate_tetrahedron:
    def __init__(self):
        self.x_1 = nucleotide()
        self.x_2 = nucleotide()
        self.x_3 = nucleotide()
        self.x_4 = nucleotide()
    
        self.e_1_2_constraint = 'NaN'
        self.e_1_3_constraint = 'NaN'
        self.e_1_4_constraint = 'NaN'
        self.e_2_3_constraint = 'NaN'
        self.e_2_4_constraint = 'NaN'
        self.e_3_4_constraint = 'NaN'
    
        self.x_1.number = -1
        self.x_2.number = -1
        self.x_3.number = -1
        self.x_4.number = -1
    
    def ntnums(self):
        return (self.x_1.number, self.x_2.number, self.x_3.number, self.x_4.number)
    
    def show(self):
        print("###")
        print("NT IDS     : %s,    %s,    %s,    %s"%(self.x_1_nt,self.x_2_nt,self.x_3_nt,self.x_4_nt))
        print("Vt IDS     : %s,    %s,    %s,    %s"%(self.x_1_n,self.x_2_n,self.x_3_n,self.x_4_n))
        print("X1 Coords  : %s"%(self.x_1))
        print("X2 Coords  : %s"%(self.x_2))
        print("X3 Coords  : %s"%(self.x_3))
        print("X4 Coords  : %s"%(self.x_4))
        print("1-2 Const  : %s"%(self.e_1_2_constraint))
        print("1-3 Const  : %s"%(self.e_1_3_constraint))
        print("1-4 Const  : %s"%(self.e_1_4_constraint))
        print("2-3 Const  : %s"%(self.e_2_3_constraint))
        print("2-4 Const  : %s"%(self.e_2_4_constraint))
        print("3-4 Const  : %s"%(self.e_3_4_constraint))
    
    def show_minimal(self):
        print self.ntnums()
        print("%s, %s, %s, %s, %s, %s, %s, %s, %s, %s"%(self.x_1.nt,self.x_2.nt,self.x_3.nt,self.x_4.nt, self.e_1_2_constraint, self.e_1_3_constraint, self.e_1_4_constraint, self.e_2_3_constraint, self.e_2_4_constraint, self.e_3_4_constraint))

def load_three_tree(tree_file, sequence):
    f = open(tree_file)
    tree_contents = f.readlines()
    tree_lines = []
    constrained_tetrahedra = []
    r = re.compile(r'(?:[^,(]|\([^)]*\))+')
    for i in tree_contents:
        split = r.findall(i)
        tree_lines.append(split)
        if(len(tree_lines) > 1):
            val = tree_lines[-1:]
            constrained_tetrahedra.append(val[0][-5:])
    rct = []
    for i in constrained_tetrahedra:
        real_i = []
        for x in i:
            x = x.replace('[',"")
            x = x.replace(']',"")
            x = x.replace('\n',"")
            x = re.split(r'[()]', x)
            for j in x:
                if(j == ' '):
                    x.remove(' ')
                elif(j == ''):
                    x.remove('')
            if(len(x) == 1):
                x[0] = x[0].replace(' ',"")
                x = x[0]
            real_i.append(x)
        rct.append(real_i)
    ThreeTree = []
    track_index = 0
    for i in rct:
        if(len(i) < 3):
            continue
        temp = candidate_tetrahedron()
        temp.x_1.number = int(i[0])
        temp.x_2.number = int(i[1])
        temp.x_3.number = int(i[2])
        temp.x_4.number = int(i[3])
        temp.x_1.nt = sequence[int(i[0])-1]
        temp.x_2.nt = sequence[int(i[1])-1]
        temp.x_3.nt = sequence[int(i[2])-1]
        temp.x_4.nt = sequence[int(i[3])-1]
        for j in i[4]:
            if(len(j)<3):
                continue
            item = j.split()
            edge_num = int(item[0]) + int(item[1])
            if(edge_num == 3):
                temp.e_1_2_constraint = item[2]
            elif(edge_num == 4):
                temp.e_1_3_constraint = item[2]
            elif(edge_num == 5 and item[0] == 1):
                temp.e_1_4_constraint = item[2]
            elif(edge_num == 5 and item[0] == 2):
                temp.e_2_3_constraint = item[2]
            elif(edge_num == 6):
                temp.e_2_4_constraint = item[2]
            elif(edge_num == 7):
                temp.e_3_4_constraint = item[2]
        if(temp.x_1.number - temp.x_2.number == 1 or temp.x_1.number-temp.x_2.number == -1):
            temp.e_1_2_constraint +="-BB-"
        if(temp.x_1.number-temp.x_3.number == 1 or temp.x_1.number-temp.x_3.number == -1):
            temp.e_1_3_constraint +="-BB-"
        if(temp.x_1.number-temp.x_4.number == 1 or temp.x_1.number-temp.x_4.number == -1):
            temp.e_1_4_constraint +="-BB-"
        if(temp.x_2.number-temp.x_3.number == 1 or temp.x_2.number-temp.x_3.number == -1):
            temp.e_2_3_constraint +="-BB-"
        if(temp.x_2.number-temp.x_4.number == 1 or temp.x_2.number-temp.x_4.number == -1):
            temp.e_2_4_constraint +="-BB-"
        if(temp.x_3.number-temp.x_4.number == 1 or temp.x_3.number-temp.x_4.number == -1):
            temp.e_3_4_constraint +="-BB-"
        ThreeTree.append(temp)
        del temp
        track_index = track_index + 1
    return ThreeTree


def database_coords_to_numpy(x):
    x=x.replace("*", "'")
    x=x.replace("^"," ")
    x = x.split()
    x = np.asarray(x)
    x.reshape(3, -1)
    x.astype(float)
    return x

def database_types_to_numpy(x):
    x=x.replace("*", "'")
    x=x.replace("^"," ")
    x = x.split()
    x = np.asarray(x)
    return x


def search_to_tetrahedron(data):
    tetra = candidate_tetrahedron()
    tetra.x_1.nt = data[0].encode('ascii')
    tetra.x_2.nt = data[1].encode('ascii')
    tetra.x_3.nt = data[2].encode('ascii')
    tetra.x_4.nt = data[3].encode('ascii')
    
    tetra.x_1.coords = database_coords_to_numpy(data[4])
    tetra.x_1.labels = database_types_to_numpy(data[5])
    tetra.x_2.coords = database_coords_to_numpy(data[6])
    tetra.x_2.labels = database_types_to_numpy(data[7])
    tetra.x_3.coords = database_coords_to_numpy(data[8])
    tetra.x_3.labels = database_types_to_numpy(data[9])
    tetra.x_4.coords = database_coords_to_numpy(data[10])
    tetra.x_4.labels = database_types_to_numpy(data[11])
    
    tetra.e_1_2_constraint = data[12]
    tetra.e_1_3_constraint = data[13]
    tetra.e_1_4_constraint = data[14]
    tetra.e_2_3_constraint = data[15]
    tetra.e_2_4_constraint = data[16]
    tetra.e_3_4_constraint = data[17]
    return tetra

#returns the candidates who match in nt types and constraints
def fetch_tetrahedron_cantidates(c, tetra, candidate_limit = 10):
    
    command  = "SELECT *, CASE WHEN nt1 IS '%s' THEN 1 ELSE 0 END "%(tetra.x_1.nt)
    command += "+ CASE WHEN nt2 IS '%s' THEN 1 ELSE 0 END "%(tetra.x_2.nt)
    command += "+ CASE WHEN nt3 IS '%s' THEN 1 ELSE 0 END "%(tetra.x_3.nt)
    command += "+ CASE WHEN nt4 IS '%s' THEN 1 ELSE 0 END "%(tetra.x_4.nt)
    
    command += "+ CASE WHEN const_1_2 IS '%s' THEN 1 ELSE 0 END "%(tetra.e_1_2_constraint)
    command += "+ CASE WHEN const_1_3 IS '%s' THEN 1 ELSE 0 END "%(tetra.e_1_3_constraint)
    command += "+ CASE WHEN const_1_4 IS '%s' THEN 1 ELSE 0 END "%(tetra.e_1_4_constraint)
    command += "+ CASE WHEN const_2_3 IS '%s' THEN 1 ELSE 0 END "%(tetra.e_2_3_constraint)
    command += "+ CASE WHEN const_2_4 IS '%s' THEN 1 ELSE 0 END "%(tetra.e_2_4_constraint)
    command += "+ CASE WHEN const_3_4 IS '%s' THEN 1 ELSE 0 END AS Matches "%(tetra.e_3_4_constraint)
    command += "FROM candidate_tetrahedra "
    command += "WHERE Matches IS NOT NULL "
    command += "ORDER BY Matches DESC "
    command += "LIMIT %s;"%(candidate_limit)
    
    c.execute(command)
    return c.fetchall()

def filter_matches(data):
    new_data = []
    for i in data:
        if(i[-1] != 10):
            data.remove(i)
            continue
        new_data.append(i)
    return new_data

  
