# Load the structures from the Zaretzki dataset and get the number of non-hydrogen atoms in each molecule
# Plot with the frequency bar chart

import os
import pandas as pd
from schrodinger.structure import StructureReader, StructureWriter
import matplotlib.pyplot as plt


import shutil
import sys
sys.path.append("../src")
from utils import renew_folder


# Load the input dataset
reader = StructureReader("../data/3A4.sdf")
dir_mae = "../data/mae"

# Get the number of non-hydrogen atoms in each molecule
num_atoms = []
molecule_rank = []
for i, structure in enumerate(reader, start=1):
    # Remove the hydrogen atoms using schrodinger,
    # index 333 5betacholestane-3_7_12_trihydroxy has hydrogen atoms not removed
    structure.deleteAtoms([atom for atom in structure.atom if atom.element == "H"])
    num_atoms.append(len(structure.atom))
    molecule_rank.append((i, len(structure.atom)))

# Step 1: Get the file names for the molecules with more than 40 atoms and less than 10 atoms
file_names = []
for f in os.listdir(dir_mae):
    idx = int(f.split("_")[1].split(".")[0])
    if idx in [
        m[0] for m in sorted(molecule_rank, key=lambda x: x[1], reverse=True)[:36]
    ]:
        file_names.append(f)
    if idx in [m[0] for m in sorted(molecule_rank, key=lambda x: x[1])[:14]]:
        file_names.append(f)

# Step 2: Get the file names for the molecules with uncommon atoms
uncommon_atoms_file_names = []
reader = StructureReader("../data/3A4.sdf")
for i, structure in enumerate(reader, start=1):
    for atom in structure.atom:
        if atom.element not in ["H", "C", "N", "O", "S", "F", "P", "Cl", "Br"]:
            uncommon_atoms_file_names.append(f"structure_{i:03d}.mae")
            break

print(f"Structures with uncommon atoms: {uncommon_atoms_file_names}")
file_names.extend(uncommon_atoms_file_names)

# Step 3: Get the file names for the molecules with peroxides
som_df = pd.read_csv('../data/3a4_all_som.csv')

peroxide_file_names = []
reader = StructureReader("../data/3A4.sdf")
for i, structure in enumerate(reader, start=1):
    for atom in structure.atom:
        if atom.index in som_df[som_df['mol_index'] == i]['som'].values:
            # check if the atom is connected to an oxygen
            # get the bonds of the atom
            som_bonds = [bond for bond in structure.bond if atom.index in [bond.atom1.index, bond.atom2.index]]
            som_oxygen_bond = None
            for som_bond in som_bonds:
                if som_bond.atom1.element == "O":
                    oxygen_atom = som_bond.atom1
                    som_oxygen_bond = som_bond
                elif som_bond.atom2.element == "O":
                    oxygen_atom = som_bond.atom2
                    som_oxygen_bond = som_bond
                else:
                    oxygen_atom = None
                if oxygen_atom is not None:
                    # get the bonds of the oxygen
                    oxygen_bonds = [bond for bond in structure.bond if oxygen_atom.index in [bond.atom1.index, bond.atom2.index]]
                    # remove the som_oxygen_bond
                    oxygen_bonds.remove(som_oxygen_bond)
                    for oxygen_bond in oxygen_bonds:
                        if oxygen_bond.atom1.element == "O" and oxygen_bond.atom2.element == "O":
                            peroxide_file_names.append(f"structure_{i:03d}.mae")
                            break
                break

peroxide_file_names = list(set(peroxide_file_names))
print(f"Structures with peroxides: {peroxide_file_names}")
file_names.extend(peroxide_file_names)

print(f"The files to be removed: {file_names}")

# Copy the remaining files to the new folder
renew_folder("../data/mae_c_clean")
renew_folder("../data/mae_n_clean")
renew_folder("../data/mae_s_clean")

# Carbon
carbon_file_names = os.listdir("../data/mae_c")
for file_name in carbon_file_names:
    if file_name not in file_names:
        shutil.copy(f"../data/mae_c/{file_name}", f"../data/mae_c_clean/{file_name}")

# Nitrogen
nitrogen_file_names = os.listdir("../data/mae_n")
for file_name in nitrogen_file_names:
    if file_name not in file_names:
        shutil.copy(f"../data/mae_n/{file_name}", f"../data/mae_n_clean/{file_name}")

# Sulfur
sulfur_file_names = os.listdir("../data/mae_s")
for file_name in sulfur_file_names:
    if file_name not in file_names:
        shutil.copy(f"../data/mae_s/{file_name}", f"../data/mae_s_clean/{file_name}")
