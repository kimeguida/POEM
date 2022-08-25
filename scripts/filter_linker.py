
#MIT License
#
#Copyright (c) 2022 LIT-University of Strasbourg
#
#Permission is hereby granted, free of charge, to any person obtaining a copy
#of this software and associated documentation files (the "Software"), to deal
#in the Software without restriction, including without limitation the rights
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following conditions:
#
#The above copyright notice and this permission notice shall be included in all
#copies or substantial portions of the Software.
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#SOFTWARE.


# Code used in
# Target-focused library design by pocket-applied computer vision and fragment deep generative linking
# Merveille Eguida, Christel Valencia, Marcel Hibert, Pascal Villa and Didier Rognan
# contacts: keguida[at]unistra.fr, rognan[at]unistra.fr





from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit import DataStructs
import re
from rdkit import RDLogger
from rdkit.Chem import AllChem
RDLogger.DisableLog('rdApp.*')

import argparse
import numpy as np
import operator




parser = argparse.ArgumentParser()
parser.add_argument('-i', '--ifile', type=str, required=True,
            help='generation_linker_descriptor.tsv')
parser.add_argument('-o', '--ofile', type=str, required=True,
            help='')
args = parser.parse_args()




def pass_filter(linker_, natoms_):
    print(Chem.rdMolDescriptors.CalcNumHeteroatoms(Chem.MolFromSmiles('*')))
    mol = Chem.MolFromSmiles(linker_)
    keep = True
    if mol is not None:
        n_aliphatic_ring = Descriptors.NumAliphaticRings(mol)
        n_aromatic_ring = Descriptors.NumAromaticRings(mol)
        n_ring = (max(n_aliphatic_ring, n_aromatic_ring))
        #print('n_aliphatic_ring', n_aliphatic_ring)
        #print('n_aromatic_ring', n_aromatic_ring)
        if n_ring < 1:
            n_heteroatoms = Chem.rdMolDescriptors.CalcNumHeteroatoms(mol) - 2 # remove attachments *
            # should be
            # no significant change in R1 selection: 141125 --> 141049
            #n_heteroatoms = Chem.rdMolDescriptors.CalcNumHeteroatoms(mol) - linker_.count('*') # remove attachments *

            if n_heteroatoms < 1:
                bonds_type = [b.GetBondType() for b in mol.GetBonds()]
                #print(bonds_type)
                if all([btype == Chem.rdchem.BondType.SINGLE for btype in bonds_type]) and natoms_ > 2:
                    keep = False

    return keep







discarded = []
with open(args.ifile, 'r') as f:
    for l in f:
        if l.startswith('gen'):
            continue
        cols = l.split('\n')[0].split('\t')
        gen = cols[0]
        gen_id = int(cols[1])
        gen_atom = int(cols[7])
        linker = cols[8]
        rot = int(cols[9])
        atom = int(cols[10])
        if gen_atom < atom: # wrong linkers including fragments
            #print(gen, '\t', linker)
            continue

  
        if rot == atom-1 and atom > 3: # purely linear linkers
            discarded.append((gen_id, gen, linker))
            #print(gen, linker)
            #continue


        if not pass_filter(linker, atom): # not purely linear but no ring or no heteroatom
                                            # branched linkers only composed of more than 3 carbon atoms
            discarded.append((gen_id, gen, linker))                    
            #print(gen, linker)
            #continue
print('discarded', len(discarded))



with open(args.ofile, 'w') as of:
    of.write('gen_id\tgen\tlinker\n')
    for gen_id, gen, linker in discarded:
        of.write('{}\t{}\t{}\n'.format(gen_id, gen, linker))

