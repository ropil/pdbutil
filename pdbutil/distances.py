#!/usr/bin/env python3
# coding=utf-8
import re
import math
'''
    pdbutil_distances, extracts a distance list from a PDB file
    Copyright (C) <2017>  <Robert Pilstål>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''

#library functions go here
class ChainError(Exception):
    pass

def matchAny(line):
    return True
    
def getAtomPositions(pdbfile, multi=False, m_type=matchAny, types=False):
    m_atom = re.compile("^ATOM")
    m_endmdl = re.compile("^ENDMDL")
    chains = {}
    curr_res = None
    for line in pdbfile:
            #stop reading if ENDMDL and one model only mode
        if not multi and m_endmdl.match(line):
            break
            
            #if we encounter an ATOM entry
        if m_atom.match(line):
                #and it is of specified atom type
            if m_type(line):
                chain = line[21]
                    #if we haven't seen this chain before, make a new chain list
                if not chain in chains:
                    chains[chain] = []
                    chains[chain].append([])
                    curr_res = line[22:26]
                    #append a new residue if this is a new residue than the last
                    #one read
                elif not line[22:26] == curr_res:
                    curr_res = line[22:26]
                    chains[chain].append([])
                    
                    #append position data for this atom to the current residue
                chains[chain][-1].append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                    #add also residue and atom types if specified
                if types:
                    chains[chain][-1][-1].append(line[17:20])
                    chains[chain][-1][-1].append(line[12:16])
    return chains

def getAtomPositionMultimodel(pdbfile, m_type=matchAny):
    m_atom = re.compile("^ATOM")
    m_model = re.compile("^MODEL")
    known_residues = {}
    residue_order = {}
    curr_res = None
    models = []
    for line in pdbfile:
            #stop reading if ENDMDL and one model only mode
        if m_model.match(line):
            models.append({})
            
            #if we encounter an ATOM entry
        if m_atom.match(line):
                #and it is of specified atom type
            if m_type(line):
                chain = line[21]
                    #we haven't seen hain before, init resname list and known residues
                if not chain in residue_order:
                    residue_order[chain] = []
                    known_residues[chain] = {}
                
                    #if we haven't seen this chain before, make a new chain list
                if not chain in models[-1]:
                    models[-1][chain] = []
                    models[-1][chain].append([])
                    curr_res = line[22:26]
                    
                    #append a new residue if this is a new residue than the last
                    #one read
                elif not line[22:26] == curr_res:
                    curr_res = line[22:26]
                    models[-1][chain].append([])
                    if not curr_res in known_residues[chain]:
                        residue_order[chain].append(curr_res)
                        known_residues[chain][curr_res] = True
                    
                    #append position data for this atom to the current residue
                models[-1][chain][-1].append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return [models, residue_order]

def gDist(a, b):
    dist = 10000.0
    r = range(3)
    for i in a:
        for j in b:
            dist = min(sum([math.pow(i[l]-j[l],2) for l in r]), dist)
    return math.sqrt(dist)

    #returns typenames as well
def gDistWithTypes(a, b):
    distances = []
    r = range(3)
    for i in a:
        for j in b:
            distances.append([sum([math.pow(i[l]-j[l],2) for l in r]), i[3], j[3], i[4], j[4]])
    distances.sort()
    distances[0][0] = math.sqrt(distances[0][0])
    return distances[0] 

def getDistance(chains, a, b, ca, cb):
    return gDist(chains[ca][a], chains[cb][b])

    #returns typenames as well
def getDistanceWithTypes(chains, a, b, ca, cb):
    return gDistWithTypes(chains[ca][a], chains[cb][b])

def printDistances(chains, ca, cb):
    for i in range(len(chains[ca])):
        for j in range(len(chains[cb])):
            print("\t".join([str(i+1), str(j+1), ca, cb, str(getDistance(chains, i, j, ca, cb))]))
            
def getDistanceMatrix(residues):
    matrix = [[None for j in range(len(residues))] for i in range(len(residues))]
    for i in range(len(residues)):
        for j in range(len(residues)):
            matrix[i][j] = gDist(residues[i], residues[j])
    return matrix
            
    #not implemented
def extractDistanceMatrixes(models, residue_order):
    chains = {}
    for chain in models[0]:
        chains[chain] = True
        
        #check so that each chain is represented in each model
    for i in range(len(models)):
        for chain in chains:
            if not chain in models[i]:
                raise ChainError
        for chain in models[i]:
            if not chain in chains:
                raise ChainError
    
    matrixes = {}
    for chain in chains:
        matrixes[chain] = []
        for i in range(len(residue_order[chain])):
            pass #Not implemented full yet
            
class DistanceObject:
    def __init__(self, ca="", a=-1, cb="", b=-1, dist=-1.0):
        self.ca = ca
        self.cb = cb
        self.a = int(a)
        self.b = int(b)
        self.dist = float(dist)
        
    def set(self, ca, a, cb, b, d):
        self.ca = ca
        self.cb = cb
        self.a = int(a)
        self.b = int(b)
        self.d = float(d)
        
    def __str__(self):
        return "\t".join([str(i) for i in [self.a + 1, self.b + 1, self.ca, self.cb, self.dist]])
    
    def __lt__(self, other):
        if isinstance(other, type(DistanceObject)):
            cmp = other.dist
        else:
            cmp = other
        if self.dist < cmp:
            return True
        return False
    
    def __gt__(self, other):
        if isinstance(other, type(DistanceObject)):
            cmp = other.dist
        else:
            cmp = other
        if self.dist > cmp:
            return True
        return False

class PdbPositionObject:
    def __init__(self):
        self.chains = {}
    def getAtomPositions(self, pdbfile, multi=False, m_type=matchAny):
        self.chains = getAtomPositions(pdbfile, multi=False, m_type=matchAny)
    def getDistance(self, a, b, ca, cb):
        return DistanceObject(ca=ca, a=a, cb=cb, b=b, dist=getDistance(self.chains, a, b, ca, cb))

#main definition for callable scripts
def main():
    import argparse
    import sys
    parser = argparse.ArgumentParser(description="Get distance matrix out of .pdb .")
    parser.add_argument("-c", nargs=1, default=["*"], metavar="string",
                        help="List of chains to read, default=* (all)")
    parser.add_argument("-m", action="store_true", default=False,
                        help="Read multiple models, default is only first")
    parser.add_argument("files", nargs="*", metavar="FILE",
                        help="Files for input")
    arguments = parser.parse_args(sys.argv[1:])
    files = arguments.files
    print_all = False
    for chain in arguments.c[0]:
        if chain == "*":
            print_all = True
        #use stdin if no supplied files
    if len(arguments.files) == 0:
        files = [sys.stdin]
    for f in files:
        infile = f
            #open file for reading if path to file specified
        if isinstance(f, type("")):
            infile = open(f, 'r')
        chains = getAtomPositions(infile, multi=arguments.m, m_type=matchAny)
        if print_all:
            for ca in chains:
                for cb in chains:
                    printDistances(chains, ca, cb)
        else:
            for ca in arguments.c[0]:
                for cb in arguments.c[0]:
                    printDistances(chains, ca, cb)
        infile.close()

    #if called from command line
if __name__ == '__main__':
    main()