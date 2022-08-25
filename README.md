

## Content

`envs/` --> conda environments <br>
`cdk8_structures/` --> target structures <br>
`aligned_fragments.tgz` available for download at [10.5281/zenodo.7023191](https://zenodo.org/deposit/7023191) --> output data after steps 1-3, input for step 4+ <br>
`scripts/` --> scripts to for library generation <br>
`output_files.tgz` available for download at [10.5281/zenodo.7023191](https://zenodo.org/deposit/7023191) --> data obtained at each step, and depedencies <br>


## Requirements
[Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html) environements with python 3.6+.


## Part 1/ Subpocket screening and fragments preparation

### 1. Alignment of fragments'supocket clouds to target pocket cloud with ProCare

sc-PDB subpockets and fragments <br>
* input structures from the [sc-PDB](http://bioinfo-pharma.u-strasbg.fr/scPDB/) database <br>
* fragmentation, cavity clouds and interactions computed with [IChem](http://bioinfo-pharma.u-strasbg.fr/labwebsite/downloads/IChem_v.5.2.9.tgz) <br>
&emsp;for more information on fragmentation, [sc-PDBFrag](http://bioinfo-pharma.u-strasbg.fr/scPDBFrag/) <br>

CDK8 pocket <br>
`<this_repo>/cdk8_structures/5hbh_cavityALL_p0-p1-p6.mol2` <br>

ProCare: [https://github.com/kimeguida/ProCare](https://github.com/kimeguida/ProCare) (state of our conda env: `envs/procare.yml`) <br>

`source <path_to_your_conda>` <br>
`conda activate procare` <br>
`cd aligned_fragments/` <br>
`python ../scripts/procare_launcher.py -s <subpocket> -t ../cdk8_structures/5hbh_cavityALL_p0-p1-p6.mol2 --transform --ligandtransform <fragment>` <br>

Outputs: aligned subpockets, fragments and `procare_scores.tsv`


### 2. Convert aligned fragments from mol2 to sdf
We used OpenEye python toolkits (state of our conda env: envs/oepython.yml): <br>
`conda deactivate` <br>
`conda activate oepython` <br>
`../scripts/convert.py <fragment>.mol2 <fragment>.sdf` <br>


### 3. Compute IChem interactions: 
`<path_to_your_ichem>/IChem ../cdk8_structures/5hbh_protein.mol2 <fragment>.mol2 > <fragment>.ifp`

<br>

## Part 2/ Selection and annotation of relevant fragments
To reproduce this step, the required data out of steps 1-to-3 were made availaible at [link](link)

### 4. Select top-scored subpockets and annotate fragments according to six predefined CDK8 areas
current directory: `aligned_fragments` containing data from steps 1-3 <br>

assignment of CDK8 areas <br>
`conda deactivate` <br>
`conda activate delinker` (state of our conda env: envs/delinker.yml)<br>
`python ../scripts/select_fragments_round1.py -f ../output_files/procare_scores.tsv -d . -p ../cdk8_structures/5hbh_protein.mol2 -c ../cdk8_structures/5hbh_cavityALL_p0-p1-p6.mol2` <br>

Outputs: <br>
`subpocket_p0_gate.list` which corresponds to GA1 <br>
`subpocket_p0_hinge.list` --> H <br>
`subpocket_p0_solv_1.list` --> SE2 <br>
`subpocket_p0_solv_2.list` --> SE1 <br>
`subpocket_p6_alphaC.list` --> AC <br>
`subpocket_p6_lys52.list` --> GA2 <br>
available in `<this_repo>/output_files/` <br>
<br>

## Part 3/ Fragments linking

### 5. Enumerate candidates for linking: pairs of fragments and atoms
`python ../scripts/linkable_fragments_round1_job.py --hinge subpocket_p0_hinge.list --gate subpocket_p0_gate.list --solv1 subpocket_p0_solv_1.list --solv2 subpocket_p0_solv_2.list --alphac subpocket_p6_alphaC.list --lys52 subpocket_p6_lys52.list` <br>

Outputs: <br>
`linkable_fragments_round1_<N>.list` with N in {0, 1, 2, 3, 4, 5, 6} available in `<this_repo>/output_files/` <br>

### 6. Linking with DeLinker
DeLinker: https://github.com/oxpig/DeLinker <br>

Feed DeLinker with connecatble candidates: <br>
`python ../scripts/delinker.py -f linkable_fragments_round1_<N>.list -p <your_path_to_DeLinker>/DeLinker/ > /dev/null` <br>
This step was distributed on computer clusters. <br>
Output: generation.smi renamed and zipped as `generation_<N>.smi.gz` with N in {0, 1, 2, 3, 4, 5, 6} available in `<this_repo>/output_files/`


### 7. Check generation success
Sometimes, no linker is generated and DeLinker might return truncated attempts.<br>
`python ../scripts/get_linker.py --file generation_<N>.smi.gz --fragsdir . --pathdelinker <your_path_to_DeLinker>/DeLinker/` <br>

Outputs: <br>
`generation_complete.smi` <br>
`generation_uncomplete.smi` <br>


### 8. Name molecules and filter with openEye Filter
SMILES were assigned IDs to keep track of the molecules infos. Filter will protonate and generate canonical SMILES different from RDKit's. <br>
`python ../scripts/index_generated_molecules.py -i generation_complete.smi -o generation_complete_indexed.smi` <br>

Output: `generation_complete_indexed.smi` <br>

`<path_to_openeye>/filter -in generation_complete_indexed.smi -out druglike_molecules.smi -fail druglike_failed_molecules.smi -filter ../output_files/filter_labo_cdk8.txt` <br>

Check that other annotations in the file did not affect how Filter processed the SMILES. In our case, we extracted the SMILES and indexes to a separate file `molecules.smi`. <br>

Output: `druglike_molecules.smi` <br>

### 9. Synthetic accessibility, descriptors, filtering:
SAscore from [https://github.com/rdkit/rdkit/blob/master/Contrib/SA_Score/sascorer.py](https://github.com/rdkit/rdkit/blob/master/Contrib/SA_Score/sascorer.py) <br>
`python ../scripts/get_sascore.py -i druglike_molecules.smi -o druglike_molecules_sascore.tsv` <br>
Output: `druglike_molecules_sascore.tsv` <br>

RDKit descriptors <br>
`python ../scripts/get_druglike_descriptors.py -i druglike_molecules.smi -o druglike_molecules_descriptors.tsv` <br>
`python ../scripts/get_linker_descriptors.py -i generation_complete_indexed.smi -o generation_linker_descriptors.tsv` <br>

Outputs: <br>
`druglike_molecules_descriptors.tsv` <br>
`generation_linker_descriptors.tsv` <br>

Clean generated linkers, remove too flexible, hydrophobic <br>
`python ../scripts/filter_linker.py -i generation_linker_descriptors.tsv -o linker_discarded.tsv` <br>

Output: `linker_discarded.tsv` <br>

### 10. Extract round 1 library
`python ../scripts/library_round1.py --descriptor druglike_molecules_descriptors.tsv --sascore druglike_molecules_sascore.tsv --discarded linker_discarded.tsv -o libr1.txt`

Output: `libr1.txt` <br>

### 11. Grow molecules in round 2 library from a selected hit
Example of hit compound 12 <br>
`python ../scripts/round2_fuse_mols.py --dl druglike_molecules.smi --gen generation_complete_indexed.smi --origin ../output_files/frag_origin.tsv --discarded linker_discarded.tsv --procare ../output_files/procare_scores.tsv` <br>

Outputs: <br>
`hit12_round2_mols.tsv` <br>
`hit12_round2_sascore_pass.tsv` <br>


### 12. Filter and generate round 2 library
RDKit descriptors <br>
`python ../scripts/get_round2_descriptors.py -i hit12_round2_mols.tsv -o hit12_round2_mols_descriptors.tsv` <br>

Output: `hit12_round2_mols_descriptors.tsv` <br>

candidates for synthesis  <br>
`python ../scripts/library_round2.py -i hit12_round2_mols_descriptors.tsv --sascore hit12_round2_sascore_pass.tsv -o libr2.txt` <br>

Output:`libr2.txt` <br>