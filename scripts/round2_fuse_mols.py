
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
from collections import defaultdict
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem import SDWriter
#from rdkit import RDLogger 
#RDLogger.DisableLog('rdApp.*')


import sascorer as sa




def xtract_procare(file_):
    scores = {}
    with open(file_, 'r') as f:
        for l in f:
            if l.startswith('Source'):
                continue
            cols = l.split()
            cfh_score = round(float(cols[4]), 4)
            id_ = cols[0]
            scores[id_] = cfh_score
    return scores




def get_sascore(rdkit_mol_):
    sascore = None
    if rdkit_mol_ is not None:
        try:
            sascore = sa.calculateScore(rdkit_mol_)
        except:
            pass
    return round(sascore, 4)




def xtract_origin(file_):
    frag_origin = {}
    with open(file_, 'r') as f:
        for l in f:
            if l.startswith('frag'):
                continue
            frag, _, _, _, origin, sp = l.split('\n')[0].split('\t')
            frag_origin[frag] = (origin, sp)
    return frag_origin



def xtract_gen_info(file_):
    gen_info = {}
    with open(file_, 'r') as f:
        for l in f:
            smi, f1, a1, f2, a2, linker, idx = l.split('\n')[0].split('\t')
            gen_info[idx] = (f1, a1, f2, a2, linker)
    return gen_info



def xtract_discarded(file_):
    discarded = []
    with open(file_, 'r') as f:
        for l in f:
            idx = l.split('\n')[0].split('\t')[0]
            discarded.append(idx)
    return set(discarded)



def xtract_druglike(file_, gen_info_, frag_origin_, discarded_):
    subpockets = {}
    dl2 = defaultdict(list) 
    with open(file_, 'r') as f:
        for l in f:
            smi, idx = l.split('\n')[0].split()
            if idx in discarded_:
                continue
            f1, a1, f2, a2, _ = gen_info_[idx]
            _, sp_1 = frag_origin_[f1]
            origin_2, sp_2 = frag_origin_[f2]
            subpockets[f1] = sp_1
            subpockets[f2] = sp_2
            dl2[(f1, a1, f2, a2)].append((smi, sp_1, sp_2, origin_2))
    return dl2, subpockets





if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--dl', type=str, required=True, 
                help='druglike_molecules.smi')
    parser.add_argument('--gen', type=str, required=True, 
                help='generation_complete_indexed.smi')
    parser.add_argument('--origin', type=str, required=True, 
                help='frag_origin.tsv')
    parser.add_argument('--discarded', type=str, required=True, 
                help='linker_discarded.tsv')
    parser.add_argument('--procare', type=str, required=True,
                help='procare_scores.tsv ')
    args = parser.parse_args()


    target_frag = ('cfh_ic_3dne_1_LL1_frag2.sdf', 'cfh_ic_2wi5_1_ZZ5_frag2.sdf')
    target_atm = (6, 1)
    #cfh_ic_3dne_1_LL1_frag2: hinge
    #cfh_ic_2wi5_1_ZZ5_frag2: lys52





    discarded = xtract_discarded(args.discarded)
    gen_info = xtract_gen_info(args.gen)
    frag_origin = xtract_origin(args.origin)
    dl2, subpockets = xtract_druglike(args.dl, gen_info, frag_origin, discarded)
    cfh_scores = xtract_procare(args.procare)



    target_sp = [subpockets[target_frag[0]], subpockets[target_frag[1]]]
    target_hinge_idx = target_sp.index('hinge')
    root_frag = target_frag[target_hinge_idx]
    excluded_atm = target_atm[target_hinge_idx]
    excluded_sp = target_sp
    del excluded_sp[target_hinge_idx]
    excluded_sp = excluded_sp[0] # should be 6

    kept_mols = []
    for key, val in dl2.items():
        f1, atm_1, f2, atm_2 = key
        if f1 == root_frag and atm_1 != excluded_atm:
            for mol in val: # duplicates with different
                smi, sp_1, sp_2, origin_2 = mol
                if sp_2 != excluded_sp:
                    kept_mols.append((smi, sp_2, f2, atm_1, origin_2))
                    # atm_1 is pyridine mapping






    ostr = "round2_smiles\tid\tseed_mol\tconnected_mol\thinge_idx\tfrag_2_sp\tfrag_2\torigin_2\tSAscore\tcfh_score\n"
    sa_ostr = "round2_smiles\tid\tseed_mol\tconnected_mol\thinge_idx\tfrag_2_sp\tfrag_2\torigin_2\tSAscore\tcfh_score\n"

    round2_mols = []
    unique_mols = {} # keep the last name and cross data later in analysis
    # common fragment 
    common_frag_mol = Chem.SDMolSupplier(root_frag)[0]
    # starting round 1 mol
    #seed_smiles = 'CCOc1ccc(cc1)OC(=O)c2cccnc2' # hit 37
    seed_smiles = 'CCOc1ccc(cc1)C(=O)Oc2cccnc2'
    mol_1 = Chem.MolFromSmiles(seed_smiles)

    # pyridine atom numbers from fragments changed in full round 1 mol
    # visualized with indexes and mapped atom indexes in full round 1 molecule to initial fragment indexes
    mapping_frag_mol = {1: 13,
                        2: 14,
                        3: 15,
                        5: 17,
                        6: 12}
                        #4 = pyridine azote, 16

    mol_id = 0
    for smi, sp_2, f2, atm_1, origin_2 in kept_mols:
        atm_1 = int(atm_1)
        mol_1_connect_idx = mapping_frag_mol[atm_1]
        mol_2 = Chem.MolFromSmiles(smi)
        if mol_2 is None:
            continue

        
        # searching index of non-pyrinine linker atom
        #   ...connected to one of the seed fragment (pyridine) atoms
        mol_2_connect_idx = None
        # get the pyridine atoms
        all_mcs_idx_2 = mol_2.GetSubstructMatches(common_frag_mol)
        if len(all_mcs_idx_2) != 1:
            continue
        # several pyridine can be in mol and we don't know the one at hinge, not traced in delinker
        # so we exclude those cases
        # when more than one pyridine are present, all are deleted by DeleteSubstructs
        mcs_idx_2 = all_mcs_idx_2[0]
        #print("mcs_idx_2", mcs_idx_2)
        for b in mol_2.GetBonds():
            begin_atom = b.GetBeginAtomIdx()
            end_atom = b.GetEndAtomIdx()
            if begin_atom in mcs_idx_2 and end_atom not in mcs_idx_2:
                #print(end_atom)
                mol_2_connect_idx = end_atom
            elif begin_atom not in mcs_idx_2 and end_atom in mcs_idx_2:
                #print(begin_atom)
                mol_2_connect_idx = begin_atom

        if mol_2_connect_idx is None:
            continue

        # remove the pyridine common fragment to create the R group
        mol_3 = AllChem.DeleteSubstructs(mol_2, common_frag_mol)

        # search new idx for linker atom in the R group
        mcs_idx_3 = mol_2.GetSubstructMatches(mol_3)
        mcs_idx_3 = mcs_idx_3[0] # take only the first solution when symetric
        #print("mcs_idx_3", mcs_idx_3)
        if mcs_idx_3 == () or mol_2_connect_idx not in mcs_idx_3:
            continue
        mol_3_connect_idx = mcs_idx_3.index(mol_2_connect_idx)



        # fuse
        combo = Chem.CombineMols(mol_1, mol_3)
        edcombo = Chem.EditableMol(combo)
        mol_3_connect_idx_combo = len(mol_1.GetAtoms()) + mol_3_connect_idx
        edcombo.AddBond(mol_1_connect_idx ,mol_3_connect_idx_combo, order=Chem.rdchem.BondType.SINGLE)
        fused_mol = edcombo.GetMol()
        #print(Chem.MolToSmiles(fused_mol), atm_1)
        try:
            Chem.SanitizeMol(fused_mol)
        except:
            print(smi, Chem.MolToSmiles(mol_2), Chem.MolToSmiles(fused_mol), atm_1)
            continue

        #round2_mols.append(mol_2)
        round2_mols.append(fused_mol)

        sascore = get_sascore(fused_mol)
        cfh_score = cfh_scores[f2.replace('cfh_', '').replace('.sdf', '_cavity4')]



        mol_id += 1
        mol_name = 'r2f_{}'.format(mol_id)


        ostr += "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(Chem.MolToSmiles(fused_mol),
                                                  mol_name,
                                                  seed_smiles,
                                                  smi,
                                                  atm_1,
                                                  sp_2,
                                                  f2,
                                                  origin_2,
                                                  sascore,
                                                  cfh_score)
        if sascore <= 3:
            sa_ostr += "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(Chem.MolToSmiles(fused_mol),
                                                  mol_name,
                                                  seed_smiles,
                                                  smi,
                                                  atm_1,
                                                  sp_2,
                                                  f2,
                                                  origin_2,
                                                  sascore,
                                                  cfh_score)
            # mol will be in one copy associated to any name
            # reconstruction of all names later
            unique_mols[Chem.MolToSmiles(fused_mol)] = (mol_name, fused_mol)


    print('round 2', len(ostr.split('\n'))-2)
    print('round 2 sascore', len(sa_ostr.split('\n'))-2)
    print('unique', len(unique_mols))




    with open('hit12_round2_mols.tsv', 'w') as of:
        of.write(ostr)

    with open('hit12_round2_sascore_pass.tsv', 'w') as of:
        of.write(sa_ostr)

