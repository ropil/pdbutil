#!/usr/bin/env python3
# coding=utf-8
import re
'''
    pdbutil_reader, reads and manipulates elements of a PDB file
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
class AtomInfoContainer:
    def __init__(self, data, bounds=[0, -1]):
        self.data = data
        self.bounds = bounds
        
    def __eq__(self, other):
        if isinstance(other, AtomInfoContainer):
            if self.data == other.data:
                return True
        else:
            if self.data == other:
                return True
        return False
    
    def __int__(self):
        return int(self.data)
        
    def read(self, line):
        self.data = line
        
    def compare(self, other):
        if self.data == other.data:
            return True
        return False
    
    def set(self, value):
        self.data = value
    
    def __str__(self):
        return str(self.data)
        
class AtomNumber(AtomInfoContainer):
    def __init__(self, data):
        AtomInfoContainer.__init__(self, data, bounds = [6,11])
        #self.bounds = [6,11]
            
    def read(self, line):
        #self.bounds = [6,11]
        self.data = int(line[6:11])
    
    def __str__(self):
        return '{0:5d}'.format(self.data)
    
class AtomType(AtomInfoContainer):
    def __init__(self, data):
        AtomInfoContainer.__init__(self, data, bounds = [12,16])
        #self.bounds = [12,16]
        
    def read(self, line):
        #self.bounds = [12,16]
        self.data = line[12:16]
    
class AtomAltloc(AtomInfoContainer):
    def __init__(self, data):
        AtomInfoContainer.__init__(self, data, bounds = [16,17])
        #self.bounds = [16,17]
        
    def read(self, line):
        #self.bounds = [16,17]
        self.data = line[16]
    
class AtomAAType(AtomInfoContainer):
    def __init__(self, data):
        AtomInfoContainer.__init__(self, data, bounds = [17,20])
        #self.bounds = [17,20]
        
    def read(self, line):
        self.data = line[17:20]

class AtomChain(AtomInfoContainer):
    def __init__(self, data):
        AtomInfoContainer.__init__(self, data, bounds = [21,22])
        #self.bounds = [21,22]

    def read(self, line):
        self.data = line[21]
        
class AtomAANumber(AtomInfoContainer):
    def __init__(self, data):
        AtomInfoContainer.__init__(self, data, bounds = [22,26])
        #self.bounds = [22,26]
        
    def read(self, line):
        self.data = int(line[22:26])
    
    def __str__(self):
        return '{0:4d}'.format(self.data)
    
class AtomICode(AtomInfoContainer):
    def __init__(self, data):
        AtomInfoContainer.__init__(self, data, bounds = [26,27])
        #self.bounds = [26,27]
        
    def read(self, line):
        #self.bounds = [26,27]
        self.data = line[26]
        
        
class AtomPosition(AtomInfoContainer):
    def __init__(self, data):
        AtomInfoContainer.__init__(self, data, bounds = [30,54])
        #self.bounds = [30,54]

    def read(self, line):
        #self.bounds = [30,54]
        self.data = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
        
    def __iter__(self):
        for pos in self.data:
            yield pos
        
    def __str__(self):
        return '{: >8.3f}{: >8.3f}{: >8.3f}'.format(self.data[0], self.data[1], self.data[2])
    
class AtomOccupancy(AtomInfoContainer):
    def __init__(self, data):
        AtomInfoContainer.__init__(self, data, bounds = [54,60])
        #self.bounds = [54,60]
        
    def read(self, line):
        #self.bounds = [54,60]
        self.data = float(line[54:60])
        
    def __str__(self):
        return '{: >6.3f}'.format(self.data)
    
class AtomTemperature(AtomInfoContainer):
    def __init__(self, data):
        AtomInfoContainer.__init__(self, data, bounds = [60,66])
        
    def read(self, line):
        #self.bounds = [60,66]
        self.data = float(line[60:66])
        
    def __str__(self):
        return '{: >6.3f}'.format(self.data)
    
class AtomElement(AtomInfoContainer):
    def __init__(self, data):
        AtomInfoContainer.__init__(self, data, bounds = [76,78])
        
    def read(self, line):
        #self.bounds = [76,78]
        self.data = line[76:78]
        
class AtomCharge(AtomInfoContainer):
    def __init__(self, data):
        AtomInfoContainer.__init__(self, data, bounds = [78,80])
        
    def read(self, line):
        #self.bounds = [78,80]
        self.data = line[78:80]
        
    def __str__(self):
        if len(self.data) != 2:
            return 2*" "
        else:
            return str(self.data)

class Atom:
    def __init__(self, number=1, type=" CA ", altloc=" ", aatype="ALA", chain="A", aanumber=1, icode=" ", position=[0.0, 0.0, 0.0], occupancy=1.0, temperature=0.0, element=" C", charge="  "):
        self.number = AtomNumber(number)
        self.type = AtomType(type)
        self.altloc = AtomAltloc(altloc)
        self.aatype = AtomAAType(aatype)
        self.chain = AtomChain(chain)
        self.aanumber = AtomAANumber(aanumber)
        self.icode = AtomICode(icode)
        self.position = AtomPosition(position)
        self.occupancy = AtomOccupancy(occupancy)
        self.temperature = AtomTemperature(temperature)
        self.element = AtomElement(element)
        self.charge = AtomCharge(charge)
        
    def __eq__(self, other):
        if isinstance(other, Atom):
            print("ATOM == ATOM: Not implemented")
        elif isinstance(other, type("")):
            if other == str(self.type):
                return True
        return False
        
    def read(self, line):
        self.number.read(line)
        self.type.read(line)
        self.altloc.read(line)
        self.aatype.read(line)
        self.chain.read(line)
        self.aanumber.read(line)
        self.icode.read(line)
        self.position.read(line)
        self.occupancy.read(line)
        self.temperature.read(line)
        self.element.read(line)
        self.charge.read(line)
        
    def setChainName(self, chainname):
        self.chain.set(chainname)
        
        #specifies B-factor for atom 
    def setTemperature(self, value):
        self.temperature.set(value)
        
        #Specifies atom number
    def setNumber(self, value):
        self.number.set(value)
        
        #Specifies residue sequence number
    def setAANumber(self, value):
        self.aanumber.set(value)
        
    def __str__(self):
        line = list(" "*80)
        line[0:6] = "ATOM  "
        for field in [self.number, self.type, self.altloc, self.aatype, self.chain, self.aanumber, self.icode, self.position, self.occupancy, self.temperature, self.element, self.charge]:
            line[field.bounds[0]:field.bounds[1]] = list(str(field))
        return "".join(line)
    
class Residue:
    def __init__(self, aatype=None, aanumber=None, chain=None):
        self.aatype = AtomAAType(aatype)
        self.aanumber = AtomAANumber(aanumber)
        self.chain = AtomChain(chain)
        self.atoms = []
        
    def __iter__(self):
        for atom in self.atoms:
            yield atom
            
    def __int__(self):
        return int(self.aanumber)
            
        #Checks only residue number and type
    def __eq__(self, other):
        if self.aanumber == other.aanumber and self.aatype == other.aatype:
            return True
        return False
    
    def __getitem__(self, index):
        return self.atoms[index]
    
    def __len__(self):
        return len(self.atoms)
    
    def append(self, other):
        if isinstance(other, Atom):
            self.atoms.append(other)
        else:
            raise TypeError("Expected Atom, got something else!")
        
        #Checks residue number, type and chain ID
    def compare(self, other):
        if self.aanumber.compare(other.aanumber) and self.aatype.compare(other.aatype) and self.chain.compare(other.chain):
            return True
        return False 
    
    def clearAtoms(self):
        self.atoms = []
        
    def setName(self, line):
            #initialize self
        self.aatype.read(line)
        self.aanumber.read(line)
        self.chain.read(line)
        
    def setChainName(self, chainname):
        self.chain.set(chainname)
        for atom in self:
            atom.setChainName(chainname)
       
        #specifies B-factor for whole residue (all participating atoms)
    def setTemperature(self, value):
        for atom in self.atoms:
            atom.setTemperature(value)
        
        #Specify the residue number 
    def setNumber(self, value):
        for atom in self:
            atom.setAANumber(value)
            
        #Retrieve the residue number
    def getNumber(self):
        return
    
    # Pop an atom of specified index
    def pop(self, index):
        return self.atoms.pop(index)
    
        #accept a reading buffer as a list of strings
    def read(self, buffer):
            #add new atoms
        for line in buffer:
            self.atoms.append(Atom())
            self.atoms[-1].read(line)
        self.setName(str(self.atoms[0]))
        #self.chain.set(str(self.atoms[0].chain))
            
    def __str__(self):
        return "\n".join([str(atom) for atom in self.atoms])
        
class Chain:
    def __init__(self, chain=None):
        self.chain = AtomChain(chain)
        self.residues = []
        
    def __iter__(self):
        for residue in self.residues:
            yield residue
            
    def __getitem__(self, index):
        return self.residues[index]
            
    def __eq__(self, other):
        if isinstance(other, Chain):
            if self.chain == other.chain:
                return True
        elif isinstance(other, AtomChain) or isinstance(other, type("")):
            if self.chain == other:
                return True
        return False
    
    def __len__(self):
        return len(self.residues)
    
    def __str__(self):
        return "\n".join([str(res) for res in self.residues])
            
    def append(self, residue):
        if isinstance(residue, Residue):
            self.residues.append(residue)
        else:
            raise TypeError
        
    def compare(self, other):
        if self.chain.compare(other.chain):
            return True
        return False
    
    # Pop a residue of specified index
    def pop(self, index):
        return self.residues.pop(index)
    
    def printTer(self):
        line = list(" "*80)
        line[0:6] = list("TER   ")
        for field in [self.residues[-1].atoms[-1].number, self.residues[-1].atoms[-1].aatype, self.residues[-1].atoms[-1].chain, self.residues[-1].atoms[-1].aanumber, self.residues[-1].atoms[-1].icode]:
            line[field.bounds[0]:field.bounds[1]] = list(str(field))
        return "".join(line)
    
    def getName(self):
        return str(self.chain)
        
    def setName(self, line):
        self.chain.read(line)
    
    def read(self, buffer):
        residueBuffers = []
        currentResidue = Residue(None)
        for line in buffer:
            nextResidue = Residue(None)
            nextResidue.setName(line)
            if not currentResidue.compare(nextResidue):
                residueBuffers.append([])
                currentResidue = Residue(None)
                currentResidue.setName(line)
            residueBuffers[-1].append(line)
            #skip empty buffer, read in all residues
        for buffer in residueBuffers:
            if len(buffer) > 0:
                self.residues.append(Residue(None))
                self.residues[-1].read(buffer)
                
    # Change identifier of chain
    def rename(self, chainname):
        self.chain.set(chainname)
        for residue in self:
            residue.setChainName(chainname)
                
    def renumber(self, start=1):
        num = start
        for residue in self:
            residue.setNumber(num)
            num += 1
            
    def sort(self, reverse=False):
        residuelist = [[int(residue), residue] for residue in self]
        residuelist.sort(reverse=reverse)
        self.residues = [res[1] for res in residuelist]



class Pdb:
    def __init__(self):
        self.chains = []
        
    def __iter__(self):
        for chain in self.chains:
            yield chain
            
    def __getitem__(self, index):
        return self.chains[index]
    
    def __len__(self):
        return len(self.chains)
                
    def __str__(self):
        output = ""
        for chain in self.chains:
            output += str(chain) + "\n" + chain.printTer() + "\n"
        return output[:-1]
    
    def read(self, infile, hetatm=False):
        m_atom = re.compile("^ATOM  ")
        m_hetatom = re.compile("^HETATM")
        chainBuffers = []
        chainOrder = []
        currentChain = Chain(None)
        
        for line in infile:
            if m_atom.match(line) or (m_hetatom.match(line) and hetatm):
                nextChain = Chain(None)
                nextChain.setName(line)
                if not currentChain.compare(nextChain):
                    chainBuffers.append([])
                    chainOrder.append(str(currentChain.chain))
                    currentChain = Chain(None)
                    currentChain.setName(line)
                    
                
                chainBuffers[-1].append(line)
            #skip empty buffer, read in all residues
        for buffer in chainBuffers:
            if len(buffer) > 0:
                self.chains.append(Chain(None))
                self.chains[-1].setName(buffer[0])
                self.chains[-1].read(buffer)
                
    def append(self, item):
        if isinstance(item, Residue):
            newchain = True
            for chain in self:
                
                if chain.chain.compare(item.chain):
                    newchain = False
                    chain.append(item)
            if newchain:
                self.append(Chain(chain=str(item.chain)))
                self.append(item)
        elif isinstance(item, Chain):
            self.chains.append(item)
        else:
            raise TypeError
        
    # Pop a chain of specified index
    def pop(self, index):
        return self.chains.pop(index)
            
    def renumber_atoms(self):
        num = 1
        for chain in self:
            for residue in chain:
                for atom in residue:
                    atom.setNumber(num)
                    num += 1
                    
    def renumber(self, start=1):
        for chain in self:
            chain.renumber(start=start)
            
    def sort(self, reverse=False):
        chainlist = [[str(chain.chain), chain] for chain in self]
        chainlist.sort(reverse=reverse)
        self.chains = [chain[1] for chain in chainlist]
           


class Models:
    def __init__(self):
        self.pdbs = []
        
    def __getitem__(self, index):
        return self.pdbs[index]
    
    def __setitem__(self, index, item):
        if not isinstance(item, Pdb) and not isinstance(item, type([])):
            raise TypeError("Item is of type " + str(type(item)))
        self.pdbs[index] = item
        
    def __iter__(self):
        for pdb in self.pdbs:
            yield pdb
    
    def __len__(self):
        return len(self.pdbs)
        
    def __str__(self):
            #set number of models in output
        output = self.printNumberOfModels()
            #add each model to output
        num = 1
        for pdb in self.pdbs:
            output += self.printModelStart(num)
            output += str(pdb) + "\n"
            output += self.printModelEnd()
            num += 1
        return output

    def assign(self, pdbs):
        self.pdbs = []
        if isinstance(pdbs, type([])):
            for pdb in pdbs:
                self.append(pdb)
        else:
            self.append(pdbs)
            
            
    def append(self, pdb):
        if isinstance(pdb, Pdb):
            self.pdbs.append(pdb)
        else:
            raise TypeError
        
    def printNumberOfModels(self, num=None):
        output = list(" "*80)
        output[:6] = "NUMMDL"
        if isinstance(num, type(None)):
            output[10:14] = "%-4d"%len(self)
        else:
            output[10:14] = "%-4d"%num
        return "".join(output) + "\n"
        
    def printModelStart(self, num):
        output = list(" "*80)
        output[0:6] = "MODEL "
        output[10:14] = '{0:4d}'.format(num)
        return "".join(output) + "\n"
        
    def printModelEnd(self):
        output = list(" "*80)
        output[0:6] = "ENDMDL"
        return "".join(output) + "\n"
    
    def renumber(self, start=1):
        for pdb in self:
            pdb.renumber(start=start)
            
    def renumber_atoms(self):
        for pdb in self:
            pdb.renumber_atoms()
    
    def read(self, infile):
            #always have one model, even if model section is lacking
        models = [[]]
        m_atom = re.compile("^ATOM  ")
        m_hetatom = re.compile("^HETATM")
        m_model_start = re.compile("^MODEL ")
        for line in infile:
            line = line.rstrip()
                #if starting a new model
            if m_model_start.match(line):
                if(len(models[-1]) > 0):
                        #append a new list of atoms
                    models.append([])
                #otherwise keep append atom entries to last model
            elif m_atom.match(line):
                models[-1].append(line)
            elif m_hetatom.match(line):
                models[-1].append(line)
            #create one pdb-object per model
        for model in models:
            pdb = Pdb()
            pdb.read(model)
            self.pdbs.append(pdb)

 

class inlineModelOperator(Models):
    def __init__(self, settings):
        self.settings = settings
        self.pdbs = 0
        
    def __len__(self):
        return self.pdbs
    
    # Auxiliary checker
    def check_done(self):
        return False
    
    # Count number of entries in multimodel .pdb
    def count(self, infile):
        num = 0
        m_model_start = re.compile("^MODEL ")
        m_atom = re.compile("^ATOM ")
        has_content = False
        for line in infile:
            #if starting a new model
            if m_model_start.match(line):
                num += 1
            elif m_atom.match(line):
                has_content = True
        # Correct for lack of MODEL entries
        if num == 0 and has_content:
            num = 1
        return num
    
    def finisher(self, pdb):
        raise IndexError("Finisher function undefined!")
    
    def operator(self, pdb):
        raise IndexError("Operator function undefined!")
        
    def publishPdb(self, model):
        self.pdbs += 1
        # Create a Pdb from model buffer
        pdb = Pdb()
        pdb.read(model)
        # Perform operation on Pdb
        self.operator(pdb)
        self.finisher(pdb)
        # Wipe model buffer
        return []
        
    def read(self, infile):
            #always have one model, even if model section is lacking
        model = []
        m_atom = re.compile("^ATOM  ")
        m_hetatom = re.compile("^HETATM")
        m_model_start = re.compile("^MODEL ")
        for line in infile:
            line = line.rstrip()
            #if starting a new model
            if m_model_start.match(line):
                if(len(model) > 0):
                    model = self.publishPdb(model)
                    # End if auxcheck returns True
                    if self.check_done():
                        break
            #otherwise keep append atom entries to last model
            elif m_atom.match(line):
                model.append(line)
            elif m_hetatom.match(line):
                model.append(line)
        if len(model) > 0:
            model = self.publishPdb(model)
            
def readFiles(files):
    models = Models()
    for f in files:
        infile = f
            #open file for reading if path to file specified
        if isinstance(f, type("")):
            infile = open(f, 'r')
        models.read(infile)
    return models
            

#main definition for callable scripts
def main():
    import argparse
    import sys
    parser = argparse.ArgumentParser(description="Read an do basic manipulations of PDB-files.")
    parser.add_argument("-renumber", nargs=1, default=[None], metavar="int",
                        help="Renumber the amino acid sequence starting from "
                        + "specified.")
    parser.add_argument("-c", nargs=1, default=["*"], metavar="CHAINS",
                        help="List of chains to keep, default=all")
    parser.add_argument("-num", nargs=1, default=["1"], metavar="INT",
                        help="Number of models to print (*=all), default=1")
    parser.add_argument("files", nargs="*", metavar="FILE",
                        help="Files for input")
    arguments = parser.parse_args(sys.argv[1:])
    files = arguments.files
    
    # Set parameters here
    chains = arguments.c[0]
    num_models = arguments.num[0]
    if not num_models == "*":
        num_models = int(num_models)
    renumber = arguments.renumber[0]
    if renumber is not None:
        renumber = int(renumber)
    
        #use stdin if no supplied files
    if len(arguments.files) == 0:
        files = [sys.stdin]
    for f in files:
        infile = f
            #open file for reading if path to file specified
        if isinstance(f, type("")):
            infile = open(f, 'r')
        models = Models()
        models.read(infile)
        if isinstance(num_models, type(int(1))):
            models.assign(models[:num_models])
        for i in range(len(models)):
            pdb = models[i]
                #append only wanted chains
            if not chains == "*":
                new_pdb = Pdb()
                for chain in pdb:
                    for c in chains:
                        if chain == c:
                            new_pdb.append(chain)
                            break
                models[i] = new_pdb
        if renumber is not None:
            models.renumber(start=renumber)
            models.renumber_atoms()
        print(models)
        infile.close()

    #if called from command line
if __name__ == '__main__':
    main()