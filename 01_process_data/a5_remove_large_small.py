# Load the structures from the Zaretzki dataset and get the number of non-hydrogen atoms in each molecule
# Plot with the frequency bar chart

import os
from schrodinger.structure import StructureReader
import matplotlib.pyplot as plt
import fnmatch


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
)  # High resolution and tight layout

# Print the molecule rank with the top 3 largest molecules
print("The molecules with 70-90 atoms:")
print(sorted(molecule_rank, key=lambda x: x[1], reverse=True)[0:3])
print("The molecules with 60-70 atoms:")
print(sorted(molecule_rank, key=lambda x: x[1], reverse=True)[3:12])
print("The molecules with 50-60 atoms:")
print(sorted(molecule_rank, key=lambda x: x[1], reverse=True)[12:18])
print("The molecules with 40-50 atoms:")
print(sorted(molecule_rank, key=lambda x: x[1], reverse=True)[18:36])
print("The molecules with 0-10 atoms:")
print(sorted(molecule_rank, key=lambda x: x[1])[:14])

# Remove the molecules with more than 50 atoms
dir_0_all_mae = "../data/zaretzki/0_all/mae"
dir_1_c_mae = "../data/zaretzki/1_c/mae"
dir_1_c_png = "../data/zaretzki/1_c/png"

# Load the structures from the Zaretzki dataset and remove the molecules with more than 40 atoms
reader = StructureReader("../data/zaretzki/0_all/3A4.sdf")

# Get the molecules which have atoms with atomic numbers larger than 40
large_atomic_number_molecules = []
for i, structure in enumerate(reader, start=1):
    high_atomic_numbers = [atom.atomic_number for atom in structure.atom if atom.atomic_number > 40]
    if high_atomic_numbers:
        large_atomic_number_molecules.append(i)
        print(f"Molecule {i} has atoms with atomic numbers larger than 40: {high_atomic_numbers}")

file_names_large_atoms = [f"structure_{i:03d}.mae" for i in large_atomic_number_molecules]


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

file_names.extend(file_names_large_atoms)

# Remove the structures including a boron
boron_file_names = []
reader = StructureReader("../data/zaretzki/0_all/3A4.sdf")
for i, structure in enumerate(reader, start=1):
    for atom in structure.atom:
        # If the atom is not a hydrogen, carbon, nitrogen, oxygen, or sulfur
        if atom.element not in ["H", "C", "N", "O", "S", "F", "P", "Cl", "Br"]:
            print(f"Structure {i} has a boron atom: {atom.element}")
            boron_file_names.append(f"structure_{i:03d}.mae")
            break

file_names.append("structure_066.mae")

# Remove the structures including a peroxide
file_names.append("structure_451.mae")

print(f"The files to be removed: {file_names}")

# Remove the files from the directories
for file_name in file_names:
    if os.path.exists(os.path.join(dir_1_c_mae, file_name)):
        os.remove(os.path.join(dir_1_c_mae, file_name))
    matching_pattern = file_name.replace(".mae", "").replace("structure_", "")
    matching_pattern = matching_pattern + "_*_som.png"
    for existing_file in os.listdir(dir_1_c_png):
        if fnmatch.fnmatch(existing_file, matching_pattern):
            os.remove(os.path.join(dir_1_c_png, existing_file))

