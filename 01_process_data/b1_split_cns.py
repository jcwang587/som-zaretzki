"""
Load the structures from the Zaretzki dataset 
and filter the structures with carbon, nitrogen and sulfur oxidation
1. The directory '../data/mae_c' will contain the structures with carbon oxidation 
including the primary, secondary, and tertiary SOM
2. The directory '../data/mae_n' will contain the structures with nitrogen oxidation 
including the primary, secondary, and tertiary SOM
3. The directory '../data/mae_s' will contain the structures with sulfur oxidation 
including the primary, secondary, and tertiary SOM
"""

import sys
import shutil

sys.path.append("../src")
from utils import renew_folder

from schrodinger.structure import StructureReader


def main():
    # Read the structures
    reader = StructureReader("../data/3A4.sdf")

    som_all = []
    mol_title_all = []

    # Iterate through the structures and get the som status
    for i, structure in enumerate(reader, start=1):
        som_list = []

        mol_title = structure.property.get("s_m_title")
        mol_title_all.append(mol_title)
        primary_som_int = structure.property.get("i_sd_PRIMARY\_SOM", "N/A")
        if primary_som_int == "N/A":
            primary_som_string = structure.property.get("s_sd_PRIMARY\_SOM", "N/A")
            if primary_som_string == "N/A":
                raise ValueError(f"No PRIMARY_SOM property found for structure {i}")
            else:
                primary_som_list = list(map(int, primary_som_string.split()))
                som_list.extend(primary_som_list)
        else:
            som_list.append(primary_som_int)

        secondary_som_int = structure.property.get("i_sd_SECONDARY\_SOM", "N/A")
        if secondary_som_int == "N/A":
            secondary_som_string = structure.property.get("s_sd_SECONDARY\_SOM", "N/A")
            if secondary_som_string == "N/A":
                pass
            else:
                secondary_som_list = list(map(int, secondary_som_string.split()))
                som_list.extend(secondary_som_list)
        else:
            som_list.append(secondary_som_int)

        tertiary_som_int = structure.property.get("i_sd_TERTIARY\_SOM", "N/A")
        if tertiary_som_int == "N/A":
            tertiary_som_string = structure.property.get("s_sd_TERTIARY\_SOM", "N/A")
            if tertiary_som_string == "N/A":
                pass
            else:
                tertiary_som_list = list(map(int, tertiary_som_string.split()))
                som_list.extend(tertiary_som_list)
        else:
            som_list.append(tertiary_som_int)

        if isinstance(som_list, int):
            som_list = [som_list]

        print(f"Structure {i}: {som_list}")
        som_all.append(som_list)

    # Re-read the structures
    reader = StructureReader("../data/3A4.sdf")

    # Iterate through the structures and print the element
    c_oxidation_idx = []
    n_oxidation_idx = []
    s_oxidation_idx = []

    for i, structure in enumerate(reader, start=1):
        som_indices = som_all[i - 1]
        for som_index in som_indices:
            element = structure.atom[som_index].element
            if element == "C":
                c_oxidation_idx.append(i)
            elif element == "N":
                n_oxidation_idx.append(i)
            elif element == "S":
                s_oxidation_idx.append(i)

    c_oxidation_idx = set(c_oxidation_idx)
    n_oxidation_idx = set(n_oxidation_idx)
    s_oxidation_idx = set(s_oxidation_idx)

    renew_folder("../data/mae_c")
    renew_folder("../data/mae_n")
    renew_folder("../data/mae_s")

    # copy the structures with carbon oxidation to the new folder based on the c_oxidation_idx
    for i in c_oxidation_idx:
        shutil.copy(
            f"../data/mae/structure_{i:03d}.mae",
            f"../data/mae_c/structure_{i:03d}.mae",
        )

    # copy the structures with nitrogen oxidation to the new folder based on the n_oxidation_idx
    for i in n_oxidation_idx:
        shutil.copy(
            f"../data/mae/structure_{i:03d}.mae",
            f"../data/mae_n/structure_{i:03d}.mae",
        )

    # copy the structures with sulfur oxidation to the new folder based on the s_oxidation_idx
    for i in s_oxidation_idx:
        shutil.copy(
            f"../data/mae/structure_{i:03d}.mae",
            f"../data/mae_s/structure_{i:03d}.mae",
        )


if __name__ == "__main__":
    main()
