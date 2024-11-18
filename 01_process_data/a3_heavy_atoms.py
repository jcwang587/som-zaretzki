# Load the structures from the Zaretzki dataset and get the number of non-hydrogen atoms in each molecule
# Plot with the frequency bar chart

from schrodinger.structure import StructureReader
import matplotlib.pyplot as plt


# Read the structures
reader = StructureReader("../data/3A4.sdf")

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
plt.xlabel("Number of Heavy Atoms", fontsize=28)
plt.ylabel("Frequency", fontsize=28)
plt.xlim(0, 90)
plt.ylim(0, 250)
plt.grid(True, which="both", linestyle="--", linewidth=0.5)

plt.xticks(fontsize=24)
plt.yticks(fontsize=24)

# Add numbers on top of each bar
for count, bin in zip(counts, bins):
    plt.text(
        bin + 5,
        count + 0.5,
        str(int(count)),
        ha="center",
        va="bottom",
        fontsize=26,
        color="black",
    )

# Save the plot
plt.savefig(
    "../data/large_molecule.png", dpi=300, bbox_inches="tight"
) 
