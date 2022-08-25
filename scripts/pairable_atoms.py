
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



def xtract_coords(file_):
    with open(file_, 'r') as f:
        sdf = f.read()
    atom_bond_block = sdf.split('V2000\n')[1].split('M  ')[0].split('\n')
    coordinates_heavy = {}
    coordinates_h = {}
    atoms_heavy = {}
    atoms_h = {}
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
                coordinates_heavy[indice] = [x, y, z]
                atoms_heavy[indice] = atm_name
            else:
                coordinates_h[indice] = [x, y, z]
                atoms_h[indice] = atm_name


        elif len(cols) == 7: # bond block
            atm1 = int(cols[0])
            atm2 = int(cols[1])
            t = int(cols[2])
            bonds[(atm1, atm2)] = t

    """for k, v in coordinates_heavy.items():
        print(k, v)"""

    return coordinates_heavy, atoms_heavy, coordinates_h, atoms_h, bonds




def hydrogen_connection(coords_h_, bonds_):
    indices_h = [i for i in coords_h_]
    #print(indices_h)
    heavy_h_connected = []
    for connection in bonds_:
        atm1 = connection[0]
        atm2 = connection[1]
        if atm1 in indices_h:
            if atm2 not in indices_h:
                heavy_h_connected.append((atm2, atm1)) # make heavy-H tuples
            """else:
                print('Error H-H bond')"""
        elif atm2 in indices_h:
            heavy_h_connected.append((atm1, atm2))
    return heavy_h_connected




def in_cone(atom_1_, exit_point_, atom_2_, angle_=(np.pi/4), h_=6):
    # angle_: 
    inside = False
    radius = np.tan(angle_) * h_
    vect_axis = np.array(exit_point_) - np.array(atom_1_)
    axis = vect_axis / np.linalg.norm(vect_axis) * h_ # norm of axis should be h_
    #print(np.linalg.norm(axis))
    proj_axis_dist = np.dot((np.array(atom_2_)-np.array(atom_1_)), axis)
    #print(proj_axis_dist)
    #if h_ >= proj_axis_dist >= 0: # means atom found in the right exit direction and cone heigth limited to 6
    if proj_axis_dist >= 0: # means atom found in the right exit direction
        vect_proj_axis = axis / h_ * proj_axis_dist
        vect_proj = (np.array(atom_2_) - np.array(atom_1_)) - vect_proj_axis

        proj_dist = np.linalg.norm(vect_proj)
        cone_radius_proj = proj_axis_dist * radius / h_

        if proj_dist <= cone_radius_proj:
            inside = True

    return inside



def calc_angle(atom_1_, exit_point_1_, atom_2_, exit_point_2_):
    vect_1 = exit_point_1_- atom_1_
    vect_2 = exit_point_2_- atom_2_
    # unit vectors
    vect_1 = vect_1/np.linalg.norm(vect_1)
    vect_2 = vect_2/np.linalg.norm(vect_2)
    

    ang = np.arccos(np.dot(vect_1, vect_2))
    return ang



def atom_pairs(coords_heavy_1_, atoms_heavy_1_, coords_h_1_, bonds_1_,
            coords_heavy_2_, atoms_heavy_2_, coords_h_2_, bonds_2_):
    heavy_h_connected_1 = hydrogen_connection(coords_h_1_, bonds_1_)
    heavy_h_connected_2 = hydrogen_connection(coords_h_2_, bonds_2_)

    pairs = {}

    for heavy_atm_1, h_atom_1 in heavy_h_connected_1:
        coord1 = np.array(coords_heavy_1_[heavy_atm_1])
        exit_point_1 =  np.array(coords_h_1_[h_atom_1])
        for heavy_atm_2, h_atom_2 in heavy_h_connected_2:
            coord2 =  np.array(coords_heavy_2_[heavy_atm_2])
            exit_point_2 =  np.array(coords_h_2_[h_atom_2])
            if in_cone(coord1, exit_point_1, coord2):
                if in_cone(coord2, exit_point_2, coord1):
                    dist = np.linalg.norm(coord2-coord1)
                    #ang = calc_angle(coord1, exit_point_1, coord2, exit_point_2)
                    pairs[((heavy_atm_1, atoms_heavy_1_[heavy_atm_1]),
                            (heavy_atm_2, atoms_heavy_2_[heavy_atm_2]))] = (dist)

    return pairs




def pairable_fragments(frag1_, frag2_, overlap_thr_=1.3):
    coord_1, _, _, _, _ = xtract_coords(frag1_)
    coord_1 = list(coord_1.values())
    coord_2, _, _, _, _ = xtract_coords(frag2_)
    coord_2 = list(coord_2.values())
    dist = []
    for c1 in coord_1:
        dist += [np.linalg.norm(np.array(c1)-np.array(c2)) for c2 in coord_2]
        n_overlap = [d < overlap_thr_ for d in dist].count(True)
    if n_overlap > 3:
        pairable = False
    else:
        pairable = True

    """
    for i, c1 in enumarate(coord_1):
            for j, c2 in enumarate(coord_2):
                dist.append((i, j, np.linalg.norm(np.array(c1)-np.array(c2))))

        pairs_overlap = [(i, j) for i, j, d in dist if d < overlap_thr_]
        
        if n_overlap > 3:
            pairable = False
        else:
            pairable = True

    """

    return pairable




if __name__ == '__main__':



    parser = argparse.ArgumentParser()
    parser.add_argument('-m1', '--molecule1', type=str, required=True,
                help='sdf file')
    parser.add_argument('-m2', '--molecule2', type=str, required=True,
                help='sdf file')
    args = parser.parse_args()


    coordinates_heavy_1, atoms_heavy_1, coordinates_h_1,atoms_h_1, bonds_1 =\
             xtract_coords(args.molecule1)
    coordinates_heavy_2, atoms_heavy_2, coordinates_h_2,atoms_h_2, bonds_2 =\
             xtract_coords(args.molecule2)

    pairs = atom_pairs(coordinates_heavy_1, atoms_heavy_1, coordinates_h_1, bonds_1,
                        coordinates_heavy_2, atoms_heavy_2, coordinates_h_2, bonds_2)

    for atm, d in pairs.items():
        mol1_file = os.path.basename(args.molecule1)
        mol2_file = os.path.basename(args.molecule2)
        print('{}\t{}-{}\t{}\t{}-{}\t{}'.format(mol1_file, *atm[0], mol2_file, *atm[1], d))

