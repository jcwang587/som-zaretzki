# The script is used for filtering the aromatic molecules.


import sys
sys.path.append('../src')
from utils import renew_folder

import os
import shutil
import pandas as pd
from schrodinger import adapter
from schrodinger.structure import StructureReader
from rdkit import Chem


renew_folder("../data/zaretzki/1_c/mae_bde_available")
renew_folder("../data/zaretzki/1_c/png_bde_available")

mae_dir = "../data/zaretzki/1_c/mae"
png_dir = "../data/zaretzki/1_c/png"

primary_som_df = pd.read_csv("../data/zaretzki/0_all/3a4_PRIMARY_som.csv")

for file in os.listdir(mae_dir):
    structure = StructureReader.read(mae_dir + "/" + file)
    mol_title = structure.property.get("s_m_title")
    primary_som = primary_som_df[primary_som_df["mol_title"] == mol_title]["PRIMARY_som"].tolist()

    rdkit_mol = adapter.to_rdkit(structure)
    aromatic_atoms = [atom.GetIdx() + 1 for atom in rdkit_mol.GetAtoms() if atom.GetIsAromatic()]

    # get the atoms that have no hydrogen neighbor
    non_hydrogen_atoms = [atom.GetIdx() + 1 for atom in rdkit_mol.GetAtoms() if atom.GetTotalNumHs(includeNeighbors=True) == 0]

    # get the atoms that have any double bond
    double_bond_atoms = [atom.GetIdx() + 1 for atom in rdkit_mol.GetAtoms() if any(bond.GetBondType() == Chem.BondType.DOUBLE for bond in atom.GetBonds())]
    
    # Check if all the primary som atoms are aromatic
    if not any(atom in non_hydrogen_atoms for atom in primary_som):
        if not any(atom in aromatic_atoms for atom in primary_som):
            if not any(atom in double_bond_atoms for atom in primary_som):
                shutil.copy(mae_dir + "/" + file, mae_dir + "_bde_available" + "/" + file)
                file_name_id_part = file.replace(".mae", "").replace("structure_", "")
                shutil.copy(png_dir + "/" + file_name_id_part + "_" + mol_title + "_som.png", 
                        png_dir + "_bde_available" + "/" + file_name_id_part + "_" + mol_title + "_som.png")
        

