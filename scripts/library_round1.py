
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





def xtract_discarded(discarded_):
    indexes = []
    with open(discarded_, 'r') as f:
        for l in f:
            if l.startswith('gen') or l == '\n':
                continue
            idx = int(l.split()[0])
            indexes.append(idx)
    return set(indexes)




def xtract_sascore(sascore_):
    sascores = {}
    unique = []
    smiles = {}
    with open(sascore_, 'r') as f:
        for l in f:
            if l.startswith('gen') or l == '\n':
                continue
            smi, idx, sa = l.split()
            # remove duplicates
            sascores[smi] = (round(float(sa), 4))
            smiles[int(idx)] = smi
    return sascores, smiles



def xtract_desc(descriptor_):
    descriptors = {}
    with open(descriptor_, 'r') as f:
        for l in f:
            if l.startswith('gen'):
                continue
            cols = l.split()
            smi = cols[0]
            tpsa = float(cols[1])
            mw = float(cols[2])
            hbd = int(cols[3])
            hba = int(cols[4])
            logp = float(cols[5])
            nrot = int(cols[6])
            natoms = int(cols[7])

            descriptors[smi] = (tpsa, mw, hbd, hba, logp, nrot, natoms)

    return descriptors




if __name__ == '__main__':
    

    import argparse


    parser = argparse.ArgumentParser()
    parser.add_argument('--descriptor', type=str, required=True,
                help='druglike_molecules_descriptors.tsv')
    parser.add_argument('--sascore', type=str, required=True,
                help='druglike_molecules_sascore.tsv')
    parser.add_argument('--discarded', type=str, required=True,
                help='linker_discarded.tsv')
    parser.add_argument('-o', '--ofile', type=str, required=False,
                help='libr1.txt')
    args = parser.parse_args()


    discarded_indexes = xtract_discarded(args.discarded)
    sascores, smiles = xtract_sascore(args.sascore)
    descriptors = xtract_desc(args.descriptor)


    selected = []
    selected_sascores = []

    n = 0
    for idx, smi in smiles.items():
        if sascores[smi] <= 3:
            selected_sascores.append(smi)
            n += 1
            if idx in discarded_indexes:
                continue
            selected.append(smi)


    #print(n)
    selected_sascores = set(selected_sascores)
    print('unique pass sascores', len(selected_sascores))

    selected = set(selected)
    print('unique selected', len(selected))


    if args.ofile is not None:
        with open(args.ofile, 'w') as of:
            for smi in selected:
                of.write('{}\n'.format(smi))