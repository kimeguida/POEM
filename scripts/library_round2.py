
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
from rdkit import Chem
from rdkit.Chem import Descriptors


parser = argparse.ArgumentParser()
parser.add_argument('-i', '--ifile', type=str, required=True,
            help='hit12_round2_mols_descriptors.tsv')
parser.add_argument('-o', '--ofile', type=str, required=True,
            help='cdk8_focused_lib_r2.txt')
parser.add_argument('--sascore', type=str, required=True,
            help='hit12_round2_sascore_pass.tsv')
args = parser.parse_args()







def xtract_pass_sascore(file_):
    sascore_selected = []
    with open(file_, 'r') as f:
        for l in f:
            if l.startswith('round2_smiles'):
                continue
            cols = l.split()
            #sascore = float(cols[-1])
            idx = cols[1]
            sascore_selected.append(idx)
    return sascore_selected



    



sascore_selected = xtract_pass_sascore(args.sascore)

selected = []
with open(args.ifile, 'r') as f:
        for l in f:
            if l.startswith('gen'):
                continue
            cols = l.split()
            gen_id = cols[0]
            smi = cols[1]
            mw = float(cols[2])
            tpsa = float(cols[3])
            logp = float(cols[4])
            hba = int(cols[5])
            hbd = int(cols[6])
            nrot = int(cols[7])
            natom = int(cols[8])
            nheteroatom = int(cols[9])
            nchiral = int(cols[10])

            if gen_id not in sascore_selected:
                continue

            if nchiral > 0:
                continue

            if nrot > 6:
                continue

            selected.append(smi)


#print('selected', len(set(selected)))
#print(selected[:10])

with open(args.ofile, 'w') as of:
    for smi in set(selected):
        of.write(f'{smi}\n')