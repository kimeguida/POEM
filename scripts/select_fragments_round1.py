
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





# select c-FH score >= 0.39 (F measure threshold for c-FPFH)
# filter, annotate fragments into six areas
# outputs:
# subpocket_p0_hinge.list           --> H (nickname in paper)
# subpocket_p0_gate.list            --> GA1
# subpocket_p6_lys52.list           --> GA2
# subpocket_p0_solv_2.list          --> SE1
# subpocket_p0_solv_1.list          --> SE2
# subpocket_p6_alphaC.list          --> AC

# the term subpocket /area used interchangeably




import argparse
import operator
import numpy as np
from collections import defaultdict
import os
import urllib
import json

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors

from sklearn.neighbors import NearestNeighbors

import matplotlib.pyplot as plt





parser = argparse.ArgumentParser()
parser.add_argument('-f', '--filescores', type=str, required=True,
            help='scoring file, e.g., procare_scores.tsv')
parser.add_argument('-d', '--dirifp', type=str, required=True,
            help='ifp directory, e.g., ../aligned_fragments')
parser.add_argument('-p', '--protein', type=str, required=True,
            help='target protein, e.g., ../cdk8_structures/5hbh_protein.mol2')
parser.add_argument('-c', '--cavity', type=str, required=True,
            help='target VolSite cavity,  e.g., ../cdk8_structures/5hbh_cavityALL_p0-p1-p6.mol2')
args = parser.parse_args()






def xtract_sdf_coords(sdf_):
    # fragments sdf coords
    with open(sdf_, 'r') as f:
        sdf = f.read()
    atom_bond_block = sdf.split('V2000\n')[1].split('M  ')[0].split('\n')
    coordinates_heavy = {}
    coordinates_h = {}
    bonds = {}
    for i, l in enumerate(atom_bond_block):
        cols = l.split()
        if len(cols) == 16: # coord block
            atm_name = cols[3]
            indice = i+1
            x = float(cols[0])
            y = float(cols[1])
            z = float(cols[2])
            if not atm_name.startswith('H'):
                coordinates_heavy[(indice, atm_name)] = np.array([x, y, z])
            else:
                coordinates_h[(indice, atm_name)] = np.array([x, y, z])

        elif len(cols) == 7: # bond block
            atm1 = int(cols[0])
            atm2 = int(cols[1])
            t = int(cols[2])
            bonds[(atm1, atm2)] = t

    return coordinates_heavy, coordinates_h, bonds




def xtract_protein_mol2_coords(mol2_):
    with open(mol2_, 'r') as f:
        mol2 = f.read()
    substr_chain = {}
    if mol2.find('@<TRIPOS>SUBSTRUCTURE') != -1:
        substr_block = mol2.split('@<TRIPOS>SUBSTRUCTURE\n')[1].split('@<TRIPOS>')[0].split('\n')
        del substr_block[-1]
        for l in substr_block:
            cols = l.split()
            res_mol2_id = cols[0]
            res = cols[1]
            res_name = res[:3]
            chain = cols[5]
            substr_chain[(res_mol2_id, res)] = chain

    coordinates = {}
    atom_block = mol2.split('@<TRIPOS>ATOM\n')[1].split('@<TRIPOS>')[0].split('\n')
    del atom_block[-1]
    for l in atom_block:
        cols = l.split()
        res_mol2_id = cols[6]
        res = cols[7]
        res_name = res[:3]
        try:
            res_num = int(res[3:])
        except:
            continue
        atm = cols[1]
        x = float(cols[2])
        y = float(cols[3])
        z = float(cols[4])
        try:
            chain = substr_chain[(res_mol2_id, res)]
            coordinates[(res_name, res_num, chain, atm)] = np.array([x, y, z])
        except KeyError:
            coordinates[(res_name, res_num, None, atm)] = np.array([x, y, z])

    return coordinates




def xtract_cavity_mol2_coords(mol2_):
    # aligned VolSite cavities
    with open(mol2_, 'r') as f:
        mol2 = f.read()
    coordinates = {}
    atom_block = mol2.split('@<TRIPOS>ATOM\n')[1].split('@<TRIPOS>')[0].split('\n')
    del atom_block[-1]
    for l in atom_block:
        cols = l.split()
        mol2_id = int(cols[0])
        ph4 = cols[1]
        atom_type = cols[5]
        x = float(cols[2])
        y = float(cols[3])
        z = float(cols[4])
        coordinates[(mol2_id, ph4, atom_type)] = np.array([x, y, z])

    return coordinates




def xtract_interactions(ifp_):
    # from IChem IFP files
    with open(ifp_, 'r') as f:
        ifp_data = f.read().split('\n')
    del ifp_data[-4:]

    ints = defaultdict(list)
    for l in ifp_data:
        cols = l.split('|')
        if len(cols) == 10:
            inter = cols[0].split()[0]
            atom = cols[1].split()[0]
            res_name, res_info = cols[3].split()
            res_num, res_chain = res_info.split('-')
            id_ = (res_name, int(res_num), res_chain, atom)
            ints[id_].append(inter)

    return ints



def xtract_scores(file_):
    # from procare_scores.tsv
    # columns: source | cfpfh | fpfh | cfh
    scores = {}
    with open(file_, 'r') as f:
        for l in  f:
            if l.startswith('Source'):
                continue
            cols = l.replace('\n', '').split('\t')
            cav = cols[0]
            cfh = float(cols[4])
            frag = 'cfh_' + cav.replace('_cavity4', '')
            scores[frag] = cfh
    return scores



def assign_subpockets(frag_coords_, ints_, protein_coords_):
    # name in code       nickname in paper
    # hinge              --> H 
    # p0_gate            --> GA1
    # p6_lys52           --> GA2
    # p0_solv_2          --> SE1
    # p0_solv_1          --> SE2
    # p6_alphaC          --> AC
    # the term subpocket /area used

    # definition of subpocket centers: (residue_name, residu_number, chain, atom_name)
    # hinge
    hinge_interactions = ['HBond_PROT', 'HBond_LIG', 'Ionic_PROT', 'Ionic_LIG']
    p0_hinge = {('ALA', 100, 'A', 'N'),
                ('ALA', 100, 'A', 'O'),
                ('ASP', 98, 'A', 'O')}

    # other areas
    subpockets = {'p0_gate': ('PHE', 97, 'A', 'CA'),
                    'p0_solv_1': ('HIS', 106, 'A', 'CE1'),
                    'p0_solv_2': ('ARG', 356, 'A', 'CZ'),
                    'p6_alphaC': ('SER', 62, 'A', 'CA'),
                    'p6_lys52': ('LYS', 52, 'A', 'NZ')}


    subpockets_cutoff = {'p0_gate': 6.0,
                        'p0_solv_1': 6.0,
                        'p0_solv_2': 6.0,
                        'p6_alphaC': 6.0,
                        'p6_lys52': 6.0}

    assigned_subpockets = []

    # check assignment to hinge first
    # search for interactions with hinge centers 
    for k, v in ints_.items():
        if k in p0_hinge:
            for interaction in v:
                if interaction in hinge_interactions:
                    assigned_subpockets.append('p0_hinge')
                    break

    # frags not assigned to hinge can be considered for other areas:
    if 'p0_hinge' not in assigned_subpockets:
        distances = []
        for subpocket, center in subpockets.items():
            subpocket_distances = []
            p_prot = protein_coords_[center]
            # get the smallest distance between fragment atoms and centers
            # is whithin cutoff, assign to that area
            for atom, coord in frag_coords_.items():
                dist = round(np.linalg.norm(p_prot-coord), 4)
                subpocket_distances.append(dist)
            distances.append((subpocket, min(subpocket_distances)))

        # keep the closest areas...
        min_subpocket_dist = min([x[1] for x in distances])
        closest_subpockets = [x for x in distances if x[1] == min_subpocket_dist]
        
        # ... whithin cutoff,
        for subpocket, dist in closest_subpockets:
            cutoff = subpockets_cutoff[subpocket]
            if dist < cutoff:
                assigned_subpockets.append(subpocket)

    #print(assigned_subpockets)
    # check for area compatibility to assign only one to a single fragment
    # according to priority rules:
    # p0_solv_1 > p0_solv_2
    # p6_lys52 > p6_alphaC
    # any other combination: ambiguous --> no assignment

    if len(assigned_subpockets) > 1:
        if list(sorted(assigned_subpockets)) == list(sorted(['p0_solv_1', 'p0_solv_2'])):
            assigned_subpockets = ['p0_solv_1']
        elif list(sorted(assigned_subpockets)) == list(sorted(['p6_alphaC', 'p6_lys52'])):
            assigned_subpockets = ['p6_lys52']
        else:
            assigned_subpockets = [None]
    elif len(assigned_subpockets) == 0:
            assigned_subpockets = [None]

    assigned_subpocket = assigned_subpockets[0] # just to get the string

    return assigned_subpocket




def coverage_alignment_on_cavity(frag_coords_, cavity_coords_, cutoff_=1.5):
    neigh = NearestNeighbors(n_neighbors=1, radius=cutoff_).fit(list(cavity_coords_.values()))
    #print(list(cavity_coords_.values()))
    #for atom, coord in frag_coords_.items():
        #print(coord)
    
    indices, distances = neigh.radius_neighbors(list(frag_coords_.values()))
    #indices, distances = neigh.kneighbors(list(frag_coords_.values()))
    n = 0
    for i, d in zip(indices, distances):
        #print(i, d)
        if not i.size == 0:
            #print('yes')
            n += 1
    coverage = round(float(n)/len(frag_coords_.values())*100, 2)
    #print(coverage)
    return coverage




def count_polar_interactions(ifp_):
    POLAR = ['HBond_PROT', 'HBond_LIG', 'Ionic_PROT', 'Ionic_LIG',
                'Aromatic_Face/Face', 'Aromatic_Edge/Face', 'Pi/Cation']
    with open(ifp_, 'r') as f:
        ifp_data = f.read().split('\n')
    del ifp_data[-4:]

    interactions = defaultdict(list)
    n = 0
    for l in ifp_data:
        cols = l.split('|')
        if len(cols) == 10:
            inter = cols[0].split()[0]
            if inter in POLAR:
                n += 1
                atom = cols[1].split()[0]
                res_name, res_info = cols[3].split()
                res_num, res_chain = res_info.split('-')
                id_ = (res_name, int(res_num), res_chain, atom)
                interactions[id_].append(inter)
    return n, interactions



def is_cofactor(sdf_):
    # Non exhaustive list from sc-PDB 
    cofactor_het = ['CO', 'ACP', 'ADP', 'AMP', 'ANP', 'ATP', 'H4B',
                    'CDP', 'CMP', 'COA', 'CAA', 'CAO', 'CTP', 'FAD',
                    'FMN','GDP', 'GMP', 'GNP', 'GTP', 'NAD', 'NAP',
                    'NAI', 'NDO', 'NDP', 'MCA', 'SCA', 'UMP']
    filename = os.path.basename(sdf_).split('_')
    het = filename[4].upper()
    #print(het)
    excluded = False
    if het in cofactor_het:
        #print('COFACTOR')
        excluded = True

    return excluded




def rule_of_tree(sdf_):
    rule_of_tree_compilance = False
    mol = Chem.SDMolSupplier(sdf_)[0]
    clogp = Descriptors.MolLogP(mol)
    mw = Descriptors.ExactMolWt(mol)
    nrot = Descriptors.NumRotatableBonds(mol)
    tpsa = Descriptors.TPSA(mol)
    nhbd = Descriptors.NumHDonors(mol)
    nhba = Descriptors.NumHAcceptors(mol)

    if (clogp <= 3 and nrot <= 3 and nhbd <= 3 and nhba <= 3 and
            mw < 300 and tpsa <= 60):
        rule_of_tree_compilance = True

    return rule_of_tree_compilance





def write_subpocket(list_sdf_, ofile_):
    with open(ofile_, 'w') as of:
        for frag in list_sdf_:
            of.write('{}\n'.format(frag))



if __name__ == '__main__':
    


    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--filescores', type=str, required=True,
                help='scoring file, e.g., procare_scores.tsv')
    parser.add_argument('-d', '--dirifp', type=str, required=True,
                help='ifp directory, e.g., ../aligned_fragments')
    parser.add_argument('-p', '--protein', type=str, required=True,
                help='target protein, e.g., ../cdk8_structures/5hbh_protein.mol2')
    parser.add_argument('-c', '--cavity', type=str, required=True,
                help='target VolSite cavity,  e.g., ../cdk8_structures/5hbh_cavityALL_p0-p1-p6.mol2')
    args = parser.parse_args()






    protein_coords = xtract_protein_mol2_coords(args.protein)
    cavity_coords = xtract_cavity_mol2_coords(args.cavity)
    cfh = xtract_scores(args.filescores)
    cfh_score_threshold = 0.39


    print("total", len(cfh))

    # counters for selection stages
    n_score = 0
    n_coverage = 0
    n_cofactor = 0
    n_ro3 = 0

    # fragments
    subpocket_p0_hinge = []
    subpocket_p0_gate = []
    subpocket_p0_solv_1 = []
    subpocket_p0_solv_2 = []
    subpocket_p6_alphaC = []
    subpocket_p6_lys52 = []


    # workflow
    for frag, score in cfh.items():
        if score >= cfh_score_threshold:
            n_score += 1
            if not is_cofactor('{}.sdf'.format(frag)):
                n_cofactor += 1
                frag_coords_heavy, frag_coords_h, _  = xtract_sdf_coords('{}.sdf'.format(frag))
                coverage_alignment = coverage_alignment_on_cavity(frag_coords_heavy, cavity_coords)
                if coverage_alignment >= 50:
                    n_coverage += 1
                    if rule_of_tree('{}.sdf'.format(frag)):
                        n_ro3 += 1
                        ints = xtract_interactions('{}/{}.ifp'.format(args.dirifp, frag))
                        pocket = assign_subpockets(frag_coords_heavy, ints, protein_coords)

                        if pocket == 'p0_hinge':
                            #print('{}.sdf'.format(frag), coverage_alignment, pocket)
                            #print('load {}.sdf'.format(frag))
                            subpocket_p0_hinge.append('{}.sdf'.format(frag))
                        if pocket == 'p0_gate':
                            subpocket_p0_gate.append('{}.sdf'.format(frag))
                        if pocket == 'p0_solv_1':
                            subpocket_p0_solv_1.append('{}.sdf'.format(frag))
                        if pocket == 'p0_solv_2':
                            subpocket_p0_solv_2.append('{}.sdf'.format(frag))
                        if pocket == 'p6_alphaC':
                            subpocket_p6_alphaC.append('{}.sdf'.format(frag))
                        if pocket == 'p6_lys52':
                            subpocket_p6_lys52.append('{}.sdf'.format(frag))



    print('pass score filter', n_score)
    print('pass cofactor filter', n_cofactor)
    print('pass coverage/buriedness filter', n_coverage)
    print('pass ro3 filter', n_ro3)



    print("subpocket_p0_hinge", len(subpocket_p0_hinge))
    print("subpocket_p0_gate", len(subpocket_p0_gate))
    print("subpocket_p0_solv_1", len(subpocket_p0_solv_1))
    print("subpocket_p0_solv_2", len(subpocket_p0_solv_2))
    print("subpocket_p6_alphaC", len(subpocket_p6_alphaC))
    print("subpocket_p6_lys52", len(subpocket_p6_lys52))



    #for frag in subpocket_p0_hinge:
    #    print('load {}'.format(frag))

    #print('\nsubpocket_p6_lys52')
    #for frag in subpocket_p6_lys52:
    #    print('load {}'.format(frag))

    #print('\nsubpocket_p0_gate')
    #for frag in subpocket_p0_gate:
    #    print('load {}'.format(frag))


    write_subpocket(subpocket_p0_hinge, 'subpocket_p0_hinge.list')
    write_subpocket(subpocket_p0_gate, 'subpocket_p0_gate.list')
    write_subpocket(subpocket_p0_solv_1, 'subpocket_p0_solv_1.list')
    write_subpocket(subpocket_p0_solv_2, 'subpocket_p0_solv_2.list')
    write_subpocket(subpocket_p6_alphaC, 'subpocket_p6_alphaC.list')
    write_subpocket(subpocket_p6_lys52, 'subpocket_p6_lys52.list')

