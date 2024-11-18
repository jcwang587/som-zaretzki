# Load the structures from the Zaretzki dataset and get the number of non-hydrogen atoms in each molecule
# Plot with the frequency bar chart

import os
import pandas as pd
from schrodinger.structure import StructureReader, StructureWriter
import matplotlib.pyplot as plt
import fnmatch

from a4_table_gen_opt import main as table_gen

import shutil
import sys
sys.path.append("../src")
from utils import renew_folder

# Run the table generation script
print("Running the SOM csv generation script...")
table_gen(dataset_path="../data/zaretzki/0_all/3A4.sdf")

# Read the structures
reader = StructureReader("../data/zaretzki/0_all/3A4.sdf")

# Get the number of non-hydrogen atoms in each molecule
num_atoms = []
molecule_rank = []
for i, structure in enumerate(reader, start=1):
    # Remove the hydrogen atoms using schrodinger,
    # index 333 5betacholestane-3_7_12_trihydroxy has hydrogen atoms not removed
    structure.deleteAtoms([atom for atom in structure.atom if atom.element == "H"])
    num_atoms.append(len(structure.atom))
    molecule_rank.append((i, len(structure.atom)))

# Plot with the distribution plot
plt.style.use("seaborn-poster")
plt.figure(figsize=(12, 8))

ax = plt.gca()
for spine in ax.spines.values():
    spine.set_linewidth(2)

counts, bins, patches = plt.hist(
    num_atoms, bins=range(0, 100, 10), color="#08312a", edgecolor="white"
)
plt.xlabel("Number of Non-Hydrogen Atoms", fontsize=18)
plt.ylabel("Frequency", fontsize=18)
plt.xlim(0, 90)
plt.grid(True, which="both", linestyle="--", linewidth=0.5)

# Add numbers on top of each bar
for count, bin in zip(counts, bins):
    plt.text(
        bin + 5,
        count + 0.5,
        str(int(count)),
        ha="center",
        va="bottom",
        fontsize=16,
        color="black",
    )

# Save the plot
plt.savefig(
    "../data/zaretzki/0_all/large_molecule.png", dpi=300, bbox_inches="tight"
)  

# Print the molecule rank with the top 3 largest molecules
print("The molecules with 0-10 atoms:")
print(sorted(molecule_rank, key=lambda x: x[1])[:14])
print("The molecules with 70-90 atoms:")
print(sorted(molecule_rank, key=lambda x: x[1], reverse=True)[0:3])
print("The molecules with 60-70 atoms:")
print(sorted(molecule_rank, key=lambda x: x[1], reverse=True)[3:12])
print("The molecules with 50-60 atoms:")
print(sorted(molecule_rank, key=lambda x: x[1], reverse=True)[12:18])
print("The molecules with 40-50 atoms:")
print(sorted(molecule_rank, key=lambda x: x[1], reverse=True)[18:36])

# Remove the molecules with more than 50 atoms
dir_0_all_mae = "../data/zaretzki/0_all/mae_all"

# Get the file names for the molecules with more than 40 atoms and less than 10 atoms
file_names = []
for f in os.listdir(dir_0_all_mae):
    idx = int(f.split("_")[1].split(".")[0])
    if idx in [
        m[0] for m in sorted(molecule_rank, key=lambda x: x[1], reverse=True)[:36]
    ]:
        file_names.append(f)
    if idx in [m[0] for m in sorted(molecule_rank, key=lambda x: x[1])[:14]]:
        file_names.append(f)

# Remove the structures including uncommon atoms
uncommon_atoms_file_names = []
reader = StructureReader("../data/zaretzki/0_all/3A4.sdf")
for i, structure in enumerate(reader, start=1):
    for atom in structure.atom:
        if atom.element not in ["H", "C", "N", "O", "S", "F", "P", "Cl", "Br"]:
            uncommon_atoms_file_names.append(f"structure_{i:03d}.mae")
            break

print(f"Structures with uncommon atoms: {uncommon_atoms_file_names}")
file_names.extend(uncommon_atoms_file_names)

# Remove the structures including SOM connected to a peroxide
som_df = pd.read_csv('../data/zaretzki/0_all/3a4_all_som.csv')

peroxide_file_names = []
reader = StructureReader("../data/zaretzki/0_all/3A4.sdf")
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