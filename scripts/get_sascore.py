
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
from rdkit import DataStructs
import sascorer as sa
import argparse


parser = argparse.ArgumentParser()
parser.add_argument('-i', '--ifile', type=str, required=True,
            help='')
parser.add_argument('-o', '--ofile', type=str, required=True,
            help='')
args = parser.parse_args()



def get_sascore(smi_):
    sascore = None
    m = Chem.MolFromSmiles(smi_)
    if m is not None:
        try:
            sascore = sa.calculateScore(m)
        except:
            pass
    return sascore



with open(args.ofile, 'w') as of:
    of.write('gen\tid\tsascore\n')
    with open(args.ifile, 'r') as f:
        for l in f:
            try:
                smi, idx = l.split()
            except ValueError:
                continue
            sascore = get_sascore(smi)
            of.write('{}\t{}\t{}\n'.format(smi, idx, sascore))

