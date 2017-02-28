#!/usr/bin/env python3
# coding=utf-8
import pdbutil.reader as pdb_reader
'''
    pdbutil_remove_atoms, removes specified atom types from PDB files
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
class Operator(pdb_reader.inlineModelOperator):
    
    def finisher(self, pdb):
        print(self.printModelStart(len(self)) + str(pdb) + '\n' +
              self.printModelEnd().rstrip('\n'))
    
    # Clear pdb of unvanted atoms
    def operator(self, pdb):
        chain_pops = []
        # For every chain
        for ci in range(len(pdb)):
            residue_pops = []
            # and residue in chain
            for ri in range(len(pdb[ci])):
                atom_pops = []
                # and atom in residue
                for ai in range(len(pdb[ci][ri])):
                    # mark for removal if not indicated as a keeper
                    if str(pdb[ci][ri][ai].type).strip() not in \
                    self.settings['keep']:
                        atom_pops.append(ai)
                # Reversely pop removed atoms
                for ai in reversed(atom_pops):
                    pdb[ci][ri].pop(ai)
                # Mark for removal if residue is empty
                if len(pdb[ci][ri]) == 0:
                    residue_pops.append(ri)
            # Reversely pop residues
            for ri in reversed(residue_pops):
                pdb[ci].pop(ri)
            # Mark chain for removal if empty
            if len(pdb[ci]) == 0:
                chain_pops.append(ci)
        # Rev pop chains
        for ci in reversed(chain_pops):
            pdb.pop(ci)
        # Renumber atoms
        pdb.renumber_atoms()

#main definition for callable scripts
def main():
    import argparse
    import sys
    parser = argparse.ArgumentParser(
        description="Removes atoms from pdb-files.")
    parser.add_argument(
        "-keep", nargs=1, default=["CA"], metavar="str[,str[,...]",
        help="Keep only specified atom types. Default=CA")
    parser.add_argument("files", nargs="*", metavar="FILE",
                        help="Files for input")
    arguments = parser.parse_args(sys.argv[1:])
    files = arguments.files
    
    #use stdin if no supplied files
    if len(arguments.files) == 0:
        files = [sys.stdin]
        
    # Set keeped atoms
    keep = {}
    for atom in arguments.keep[0].split(','):
        keep[atom] = 1
    modelreader = Operator({'keep' : keep})
    num = 0
    
    # Count the models
    for f in files:
        infile = f
            #open file for reading if path to file specified
        if isinstance(f, type("")):
            infile = open(f, 'r')
            num += modelreader.count(infile)
            infile.close()
        else:
            num += modelreader.count(infile)
            infile.seek(0)
    
    print(modelreader.printNumberOfModels(num).rstrip('\n'))
    
    # Print the models
    for f in files:
        infile = f
            #open file for reading if path to file specified
        if isinstance(f, type("")):
            infile = open(f, 'r')
            modelreader.read(infile)
            infile.close()
        else:
            modelreader.read(infile)

    #if called from command line
if __name__ == '__main__':
    main()