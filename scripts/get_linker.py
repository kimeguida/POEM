
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



import argparse
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
import gzip
from rdkit import RDLogger 
RDLogger.DisableLog('rdApp.*')




parser = argparse.ArgumentParser()
parser.add_argument('--file', type=str, required=True,
            help='generation file, e.g., generation_0.smi.gz')
parser.add_argument('--fragsdir', type=str, required=True,
                help='fragments directory')
parser.add_argument('--pathdelinker', type=str, required=True,
                help='your path to DeLinker repo')
args = parser.parse_args()



import sys
sys.path.append('{}/analysis/'.format(args.pathdelinker))
#sys.path.append('/projects/PointCloud/DeLinker/')

from frag_utils import get_linker



def xtract_gen_info(zip_file_, fragsdir):
    infos = []
    with gzip.open(zip_file_, 'rb') as f:
        for l in f:
            #print(l)
            l = l.decode('utf-8')
            cols = l.split('\t')
            gen_smi = cols[0]
            frag_1 = cols[1]
            idx_1 = int(cols[2])
            frag_2 = cols[3]
            idx_2 = int(cols[4])
            frag_1_sdf = Chem.SDMolSupplier('{}/{}'.format(fragsdir, frag_1))[0]
            frag_2_sdf = Chem.SDMolSupplier('{}/{}'.format(fragsdir, frag_2))[0]
            combo_no_exit = Chem.CombineMols(frag_1_sdf, frag_2_sdf)
            combo = Chem.CombineMols(combo_no_exit, Chem.MolFromSmiles("*.*"))

            combo_idx_1 = idx_1 -1 # rdkit indices starts with 0
            combo_idx_2 = len(frag_1_sdf.GetAtoms()) + idx_2 -1
            #print(combo_idx_1, combo_idx_2)

            edcombo = Chem.EditableMol(combo)
            num_heavy_atoms = combo.GetNumHeavyAtoms()
            edcombo.AddBond(num_heavy_atoms, combo_idx_1,
                                order=Chem.rdchem.BondType.SINGLE)
            edcombo.AddBond(num_heavy_atoms+1, combo_idx_2,
                                order=Chem.rdchem.BondType.SINGLE)
            mol_to_link = edcombo.GetMol()
            Chem.SanitizeMol(mol_to_link)


            frags_mol = mol_to_link
            frags_smi = Chem.MolToSmiles(frags_mol)
            gen_mol = Chem.MolFromSmiles(gen_smi)
            du = Chem.MolFromSmiles('*')
            clean_frags = Chem.RemoveHs(AllChem.ReplaceSubstructs(frags_mol,du,Chem.MolFromSmiles('[H]'),True)[0])
            #print(Chem.MolToSmiles(clean_frags), frags_smi )
            infos.append((gen_mol, clean_frags, frags_smi, frag_1, idx_1, frag_2, idx_2, gen_smi))

    return infos

            





uncomplete = []
complete = []

of_complete = open('generation_complete.smi', 'a')
of_uncomplete = open('generation_uncomplete.smi', 'a')

infos = xtract_gen_info(args.file, args.fragsdir)
for gen_mol, clean_frags, frags_smi, frag_1, idx_1, frag_2, idx_2, gen_smi in infos:
    linker = get_linker(gen_mol, clean_frags, frags_smi)
    if linker == '':
        #uncomplete.append((gen_smi, frag_1, idx_1, frag_2, idx_2))
        of_uncomplete.write('{}\t{}\t{}\t{}\t{}\n'.format(gen_smi, frag_1, idx_1, frag_2, idx_2))
    else:
        #complete.append((gen_smi, frag_1, idx_1, frag_2, idx_2, linker))
        of_complete.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(gen_smi, frag_1, idx_1, frag_2, idx_2, linker))


of_complete.close()
of_uncomplete.close()

