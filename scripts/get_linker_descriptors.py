
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
import argparse
from rdkit import RDLogger 
RDLogger.DisableLog('rdApp.*')


parser = argparse.ArgumentParser()
parser.add_argument('-i', '--ifile', type=str, required=True,
            help='generation_complete_indexed.smi')
parser.add_argument('-o', '--ofile', type=str, required=True,
            help='')
args = parser.parse_args()



with open(args.ofile, 'w') as of:
    of.write('gen\tid\tf1\ta1\tf2\ta2\tgen_nrot\tgen_numatoms\tlink\tlink_nrot\tlink_numatoms\n')
    with open(args.ifile, 'r') as f:
        for l in f:
            cols = l.split()
            #print(cols)
            gen_smi, f1, a1, f2, a2, link_smi, idx = cols
            m_link = Chem.MolFromSmiles(link_smi)
            if m_link is None:
                continue
            m_gen = Chem.MolFromSmiles(gen_smi)
            if m_gen is None:
                continue
            nrot_link = Descriptors.NumRotatableBonds(m_link)
            #num_link = m_link.GetNumAtoms() - 2 # remove attachments
            num_link = m_link.GetNumAtoms() - link_smi.count('*') # remove attachments, some linkers are multicomponent *

            nrot_gen = Descriptors.NumRotatableBonds(m_gen)
            num_gen = m_gen.GetNumAtoms()

            of.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                gen_smi, idx, f1, a1, f2, a2, nrot_gen, num_gen,
                link_smi, nrot_link, num_link))

