# This script is the first step processing the Zaretzki dataset.
# The processed dataset is stored in the '/data/zaretzki/mae' folder.
# The script may take a while to run (around 30 minutes).


import sys

sys.path.append("../src")
from utils import renew_folder

from schrodinger.structure import StructureReader, StructureWriter
from schrodinger.structutils.build import add_hydrogens
from schrodinger.forcefield.minimizer import minimize_structure


# Read the structures
reader = StructureReader("../data/zaretzki/0_all/3A4.sdf")

# Count the number of structures
structure_count = sum(1 for _ in reader)

# Print the number of structures
print(f"The number of structures is: {structure_count}")

# Create a new folder to store the structures
renew_folder("../data/zaretzki/0_all/mae")

# Reset the reader to the beginning of the file
reader = StructureReader("../data/zaretzki/0_all/3A4.sdf")

# Export all the structures to multiple MAE files
for i, structure in enumerate(reader, start=1):
    print(f"Processing structure {i}...")
    add_hydrogens(structure)

    minimize_structure(structure)

    mae_path = f"../data/zaretzki/0_all/mae/structure_{i:03d}.mae"
    with StructureWriter(mae_path) as writer:
        writer.append(structure)