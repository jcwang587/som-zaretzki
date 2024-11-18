# Load the structures from the Zaretzki dataset and filter the structures with nitrogen and sulfur oxidation
# The directory '/data/zaretzki/mae_n' will contain the structures with nitrogen oxidation
# The directory '/data/zaretzki/mae_s' will contain the structures with sulfur oxidation


import sys

sys.path.append("../src")
from utils import renew_folder

import shutil
from schrodinger.structure import StructureReader

import matplotlib.pyplot as plt
from matplotlib_venn import venn3, venn3_circles

# Read the structures
reader = StructureReader("../data/zaretzki/0_all/3A4.sdf")

primary_som_all = []
mol_title_all = []

# Iterate through the structures and get the PRIMARY_SOM property
for i, structure in enumerate(reader, start=1):
    mol_title = structure.property.get("s_m_title")
    mol_title_all.append(mol_title)
    primary_som_int = structure.property.get("i_sd_PRIMARY\_SOM", "N/A")
    if primary_som_int == "N/A":
        primary_som_string = structure.property.get("s_sd_PRIMARY\_SOM", "N/A")
        if primary_som_string == "N/A":
            print("Warning")
        else:
            primary_som_list = list(map(int, primary_som_string.split()))
            primary_som_all.append(primary_som_list)
    else:
        primary_som_all.append([primary_som_int])


# Re-read the structures
reader = StructureReader("../data/zaretzki/0_all/3A4.sdf")

# Iterate through the structures and print the element
n_oxidation_idx = []
s_oxidation_idx = []
c_oxidation_idx = []
o_oxidation_idx = []
p_oxidation_idx = []
br_oxidation_idx = []

for i, structure in enumerate(reader, start=1):
    som_indices = primary_som_all[i - 1]
    for som_index in som_indices:
        element = structure.atom[som_index].element
        if element == "N":
            n_oxidation_idx.append(i)
        elif element == "S":
            s_oxidation_idx.append(i)
        elif element == "C":
            c_oxidation_idx.append(i)
        elif element == "O":
            o_oxidation_idx.append(i)
        elif element == "P":
            p_oxidation_idx.append(i)
        elif element == "Br":
            br_oxidation_idx.append(i)


n_oxidation_idx = set(n_oxidation_idx)
s_oxidation_idx = set(s_oxidation_idx)
c_oxidation_idx = set(c_oxidation_idx)
o_oxidation_idx = set(o_oxidation_idx)
p_oxidation_idx = set(p_oxidation_idx)
br_oxidation_idx = set(br_oxidation_idx)

print(o_oxidation_idx)
print(p_oxidation_idx)
print(br_oxidation_idx)


# Create a folder for nitrogen oxidation molecules
renew_folder("../data/zaretzki/1_c/mae")
renew_folder("../data/zaretzki/1_c/png")
# copy the structures with carbon oxidation to the new folder based on the c_oxidation_idx
for i in c_oxidation_idx:
    shutil.copy(
        f"../data/zaretzki/0_all/mae/structure_{i:03d}.mae",
        f"../data/zaretzki/1_c/mae/structure_{i:03d}.mae",
    )

# copy the png files to the new folder based on the c_oxidation_idx
for i in c_oxidation_idx:
    shutil.copy(
        f"../data/zaretzki/0_all/png/{i:03d}_{mol_title_all[i-1]}_som.png",
        f"../data/zaretzki/1_c/png/{i:03d}_{mol_title_all[i-1]}_som.png",
    )

# generate venn diagram for the three sets
plt.figure(figsize=(10, 8))
v = venn3(
    [n_oxidation_idx, s_oxidation_idx, c_oxidation_idx],
    ("Nitrogen", "Sulfur", "Carbon"),
)

print(len(n_oxidation_idx), len(s_oxidation_idx), len(c_oxidation_idx))
c = venn3_circles(
    subsets=[n_oxidation_idx, s_oxidation_idx, c_oxidation_idx],
    linestyle="dashed",
    linewidth=2,
)

# increase the font size of the labels and the numbers
for text in v.set_labels:
    text.set_fontsize(24)
    text.set_fontweight("bold")

# Increase the font size of the numbers inside the circles
for text in v.subset_labels:
    if text is not None:
        text.set_fontsize(24)
        text.set_fontweight("bold")

# save the venn diagram to the current directory
plt.savefig("../data/zaretzki/0_all/venn_diagram.png", dpi=300, bbox_inches="tight")
