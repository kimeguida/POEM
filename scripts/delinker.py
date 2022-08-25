
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




# customized from https://github.com/oxpig/DeLinker/blob/master/examples/DeLinker_fragment_linking.ipynb

import argparse
import sys
import os

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem.Draw import MolDrawing, DrawingOptions
from rdkit.Chem import MolStandardize
from rdkit import RDLogger 
RDLogger.DisableLog('rdApp.*')

import numpy as np

from itertools import product
import joblib
from joblib import Parallel, delayed
import re
from collections import defaultdict


parser = argparse.ArgumentParser()
parser.add_argument('-f', '--file', type=str, required=True,
            help='file of linking info: pairable atoms')
parser.add_argument('-p', '--pathdelinker', type=str, required=True,
            help='root path to delinker files')
args = parser.parse_args()


# args['pathdelinker'] = '<path_to_delinker>/DeLinker/'
# DeLinker config
sys.path.append(args.pathdelinker)
sys.path.append("{}analysis/".format(args.pathdelinker))
sys.path.append("{}examples/".format(args.pathdelinker))

# How many cores for multiprocessing
n_cores = 4
# Whether to use GPU for generating molecules with DeLinker
use_gpu = True

if not use_gpu:
    os.environ['CUDA_VISIBLE_DEVICES'] = '-1'


from DeLinker_test import DenseGGNNChemModel
import frag_utils
import rdkit_conf_parallel
from data.prepare_data import read_file, preprocess
import example_utils






def xtract_linking_info(file_):
    pairs = []
    with open(file_) as f:
        for l in f:
            cols = l.replace('\n', '').split('\t')
            frag_1 = cols[0]
            indice_1, atom_1 = cols[1].split('-')
            frag_2 = cols[2]
            indice_2, atom_2 = cols[3].split('-')
            data = (frag_1, int(indice_1), frag_2, int(indice_2))
            if data not in pairs: # duplicates in file
                pairs.append(data)
    return pairs


    

def run_delinker(sdf_1_, indice_1_, sdf_2_, indice_2_, output_prefix_):

    # read prepare input molecules
    frag_1_file = sdf_1_ #'cfh_ic_1hzz_2_NAP_frag2.sdf'
    frag_2_file = sdf_2_ #'cfh_ic_1oi9_1_N20_frag2.sdf'
    idx_1 = indice_1_ #1
    idx_2 = indice_2_ #5
    frag_1_sdf = Chem.SDMolSupplier(frag_1_file)[0]
    frag_2_sdf = Chem.SDMolSupplier(frag_2_file)[0]


    combo_no_exit = Chem.CombineMols(frag_1_sdf, frag_2_sdf)
    combo = Chem.CombineMols(combo_no_exit, Chem.MolFromSmiles("*.*"))
    #combo_2d = Chem.Mol(combo)
    #tmp = AllChem.Compute2DCoords(combo_2d)
    #example_utils.mol_with_atom_index(combo_2d)

    combo_idx_1 = idx_1 -1 # rdkit indices starts with 0
    combo_idx_2 = len(frag_1_sdf.GetAtoms()) + idx_2 -1
    #print(combo_idx_1, combo_idx_2)

    edcombo = Chem.EditableMol(combo)
    num_heavy_atoms = combo.GetNumHeavyAtoms()
    edcombo.AddBond(num_heavy_atoms, combo_idx_1,
                        order=Chem.rdchem.BondType.SINGLE)
    edcombo.AddBond(num_heavy_atoms+1, combo_idx_2,
                        order=Chem.rdchem.BondType.SINGLE)
    #editedcombo = edcombo.GetMol()
    #AllChem.Compute2DCoords(editedcombo)
    #Chem.SanitizeMol(editedcombo)
    #editedcombo

    mol_to_link = edcombo.GetMol()
    Chem.SanitizeMol(mol_to_link)

    # Convert exit vectors to carbons for conformer generation
    du = Chem.MolFromSmiles('*')
    mol_to_link_carbon = AllChem.ReplaceSubstructs(mol_to_link,
                        du, Chem.MolFromSmiles('C'),True)[0]
    Chem.SanitizeMol(mol_to_link_carbon)

    # Generate conformer
    mol_to_link_carbon = Chem.AddHs(mol_to_link_carbon)
    AllChem.ConstrainedEmbed(mol_to_link_carbon, combo_no_exit,
                        randomseed=42)
    mol_to_link_carbon = Chem.RemoveHs(mol_to_link_carbon)

    # Add this conformer to the two unlinked fragments
    conf = mol_to_link.GetConformer()
    ref_conf = mol_to_link_carbon.GetConformer()
    for i in range(mol_to_link_carbon.GetNumAtoms()):
        pos = list(ref_conf.GetAtomPosition(i))
        conf.SetAtomPosition(i, pos)
    conf.SetId(0)
    _ = mol_to_link.AddConformer(conf)


    # Get distance and angle between fragments
    dist, ang = frag_utils.compute_distance_and_angle(mol_to_link, "",
                        Chem.MolToSmiles(mol_to_link))
    Chem.MolToSmiles(mol_to_link), dist, ang
    # Write data to file
    data_path = "./{}.txt".format(output_prefix_)
    with open(data_path, 'w') as f:
        f.write("%s %s %s" % (Chem.MolToSmiles(mol_to_link), dist, ang))
    raw_data = read_file(data_path)
    preprocess(raw_data, "zinc", output_prefix_, True) # check the value of 'test' parameter False/True
    # preprocess(raw_data, dataset, name, test_mode)
    # if test_mode is true, full molecules = fragments


    # Arguments for DeLinker
    delinker_args = defaultdict(None)
    delinker_args['--dataset'] = 'zinc'
    delinker_args['--config'] = ('{"generation": true, ')+\
                       (' "batch_size": 1, ')+\
                       ('"number_of_generation_per_valid": 10, ')+\
                       ('"min_atoms": 1, "max_atoms": 6, ')+\
                       ('"train_file": "molecules_{}.json", '.format(output_prefix_))+\
                       ('"valid_file": "molecules_{}.json", '.format(output_prefix_))+\
                       ('"output_name": "{}_generation.smi"'.format(output_prefix_))+\
                       '}'

    delinker_args['--freeze-graph-model'] = False
    delinker_args['--restore'] = '{}models/pretrained_DeLinker_model.pickle'.format(
                                                                    args.pathdelinker)

    # Setup model and generate molecules
    model = DenseGGNNChemModel(delinker_args)
    model.train()
    # Free up some memory
    model = ''



def process_mol(generation_file_, info_, ofile_):
    with open(generation_file_, 'r') as f:
        smiles_unique = []
        for l in f:
            cols = l.replace('\n', '').split()
            mol = cols[2]
            if mol not in smiles_unique:
                smiles_unique.append(mol)

    with open(ofile_, 'a') as of:
        for m in smiles_unique:
            of.write('{}\t{}\t{}\t{}\t{}\n'.format(m, *info_))



if __name__ == '__main__':
    

    linking_info = xtract_linking_info(args.file)
    for info in linking_info:
        try:
            prefix = '{}_{}_{}_{}'.format(*info).replace('.sdf', '')
            run_delinker(*info, prefix)
            os.system('rm *params_zinc.json *_generated_smiles_zinc')
            process_mol('{}_generation.smi'.format(prefix), info, 'generation.smi')
            os.system('rm {0}.txt molecules_{0}.json {0}_generation.smi'.format(prefix))
        except:
            os.system('rm *params_zinc.json *_generated_smiles_zinc')
            prefix = '{}_{}_{}_{}'.format(*info).replace('.sdf', '')
            os.system('rm {0}.txt molecules_{0}.json'.format(prefix))

