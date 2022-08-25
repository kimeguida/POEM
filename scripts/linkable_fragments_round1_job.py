
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
import numpy as np
import os
import pairable_atoms as pa
import random





def read_subpocket_file(file_):
    with open(file_, 'r') as f:
        frags = f.read().split('\n')
    del frags[-1]
    return frags




def pair_fragments(frag_1_, pocket_1_, frag_2_, pocket_2_):
    coordinates_heavy_1, atoms_heavy_1, coordinates_h_1,atoms_h_1,\
            bonds_1 = pa.xtract_coords(frag_1_)
    coordinates_heavy_2, atoms_heavy_2, coordinates_h_2,atoms_h_2,\
            bonds_2 = pa.xtract_coords(frag_2_)

    connections = pa.atom_pairs(coordinates_heavy_1, atoms_heavy_1, coordinates_h_1,
                bonds_1, coordinates_heavy_2, atoms_heavy_2, coordinates_h_2,
                bonds_2)
    pairs = []
    mol1_file = os.path.basename(frag_1_)
    mol2_file = os.path.basename(frag_2_)
    for atm, d in connections.items():
        pairs.append('{}\t{}-{}\t{}\t{}-{}\t{}\t{}\t{}'.format(mol1_file, *atm[0],
                        mol2_file, *atm[1], d, pocket_1_, pocket_2_))

    return pairs





if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--hinge', type=str, required=True,
                help='file: subpocket list of fragments')
    parser.add_argument('--gate', type=str, required=True,
                help='file: subpocket list of fragments')
    parser.add_argument('--solv1', type=str, required=True,
                help='subpocket list of fragments')
    parser.add_argument('--solv2', type=str, required=True,
                help='subpocket list of fragments')
    parser.add_argument('--alphac', type=str, required=True,
                help='subpocket list of fragments')
    parser.add_argument('--lys52', type=str, required=True,
                help='subpocket list of fragments')
    args = parser.parse_args()



    #linkable subpockets in this version:
    # hinge versus solv1, lys52, gate



    # sized inputs for jobs on computer clusters
    time_per_pair = 15 # seconds
    time_total = 36 * 3600
    max_per_job = int(time_total / time_per_pair)
    print("max per job", max_per_job)



    hinge_frags = read_subpocket_file(args.hinge)
    solv1_frags = read_subpocket_file(args.solv1)
    lys52_frags = read_subpocket_file(args.lys52)
    gate_frags = read_subpocket_file(args.gate)
    #solv2_frags = read_subpocket_file(args.solv2)
    #alphac_frags = read_subpocket_file(args.alphac)
    
    # sample only 100 from the 936 fragments in lys52 (a.k.a. GA2)
    random.seed(42)
    lys52_frags_sample = random.sample(lys52_frags, 100)
    #print('lys52 sampled', len(lys52_frags_sample))

    
    # connect fragments when their positioning allow
    # the current version is flexible and might allow a few overlapping pairs
    # that can be filtered out in later stages
    # to be more stringent on the pairing and reduce the job size
    # the commented code "if pa.pairable_fragments(frag1, frag2, 2)..."
    # can be used instead


    frag_connection = []
    #pairs 

    for frag1 in hinge_frags:
        for frag2 in gate_frags:
            #if pa.pairable_fragments(frag1, frag2, 2):
            #    pairs = pair_fragments(frag1, 'hinge', frag2, 'gate')
            #    frag_connection += pairs
            pairs = pair_fragments(frag1, 'hinge', frag2, 'gate')
            frag_connection += pairs
            


        for frag2 in solv1_frags:
            #if pa.pairable_fragments(frag1, frag2, 2):
            #    pairs = pair_fragments(frag1, 'hinge', frag2, 'solv1')
            #    frag_connection += pairs
            pairs = pair_fragments(frag1, 'hinge', frag2, 'solv1')
            frag_connection += pairs

        for frag2 in lys52_frags_sample:
            #if pa.pairable_fragments(frag1, frag2, 2):
            #    pairs = pair_fragments(frag1, 'hinge', frag2, 'lys52')
            #    frag_connection += pairs
            pairs = pair_fragments(frag1, 'hinge', frag2, 'lys52')
            frag_connection += pairs


    # check content
    # format: tab-separated
    #         frag1 | frag1_atom_index-atom_symbol |
    #         frag2 | frag2_atom_index-atom_symbol | distance |
    #         area1 | area2
    # if allowed in pairable_atoms script, (distance, exit angle)

    print(len(frag_connection))
    for e in frag_connection[-10:]:
        print(e)


    # output files
    ofile = "linkable_fragments_round1"


    idx = 0
    n = 0
    for i, info in enumerate(frag_connection):
        idx = i//max_per_job
        #print(idx)
        with open('{}_{}.list'.format(ofile, idx), 'a') as of:
            of.write('{}\n'.format(info))
