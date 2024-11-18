# This script is used to generate a table of the primary SOM for the 3A4 dataset
# The csv data includes the mol_title, mol_index, primary_som, and element
# -----------------------------------------------------------
# mol_title: the title of the molecule from the Zaretzki dataset
# mol_index: the index of the molecule in the Zaretzki dataset
# primary_som: the index of the primary SOM in the molecule (could be multiple)
# element: the element of the atom in the primary SOM


from schrodinger.structure import StructureReader
import pandas as pd


def som_table_gen(som_level="PRIMARY"):
    # Read the structures
    reader = StructureReader("../data/zaretzki/0_all/3A4.sdf")

    som_all = []
    mol_title_all = []
    # Iterate through the structures and get the PRIMARY_SOM property
    for i, structure in enumerate(reader, start=1):
        mol_title = structure.property.get("s_m_title")
        mol_title_all.append(mol_title)
        som_int = structure.property.get(f"i_sd_{som_level}\_SOM", "N/A")
        if som_int == "N/A":
            som_string = structure.property.get(f"s_sd_{som_level}\_SOM", "N/A")
            if som_string == "N/A":
                som_all.append([])
            else:
                som_list = list(map(int, som_string.split()))
                som_all.append(som_list)
        else:
            som_all.append([som_int])


    # Re-read the structures
    reader = StructureReader("../data/zaretzki/0_all/3A4.sdf")

    # Get the length of primary_som_all including the length of the list in each element
    som_all_len = sum([len(som) for som in som_all])

    # initialize the pandas dataframe
    df = pd.DataFrame(columns=["mol_title", "mol_index", f"{som_level}_som", "element"])

    # Add the mol_title and primary_som to the dataframe
    for idx, (mol_title, som, structure) in enumerate(
        zip(mol_title_all, som_all, reader), start=1
    ):
        for som in som:
            element = structure.atom[som].element
            df = df.append(
                {
                    "mol_title": mol_title,
                    "mol_index": idx,
                    f"{som_level}_som": som,
                    "element": element,
                },
                ignore_index=True,
            )

    # export the dataframe to a csv file
    df.to_csv(f"../data/zaretzki/0_all/3a4_{som_level}_som.csv", index=False)


if __name__ == "__main__":  
    som_table_gen(som_level="PRIMARY")
    som_table_gen(som_level="SECONDARY")
    som_table_gen(som_level="TERTIARY")

    # concatenate the three tables
    primary_som = pd.read_csv('../data/zaretzki/0_all/3a4_PRIMARY_som.csv')
    secondary_som = pd.read_csv('../data/zaretzki/0_all/3a4_SECONDARY_som.csv')
    tertiary_som = pd.read_csv('../data/zaretzki/0_all/3a4_TERTIARY_som.csv')
    
    # change the column names
    primary_som.rename(columns={'PRIMARY_som': 'som'}, inplace=True)
    secondary_som.rename(columns={'SECONDARY_som': 'som'}, inplace=True)
    tertiary_som.rename(columns={'TERTIARY_som': 'som'}, inplace=True)

    # concatenate the three tables
    all_som = pd.concat([primary_som, secondary_som, tertiary_som], ignore_index=True)

    # sort the table by mol_index and then som
    all_som = all_som.sort_values(by=['mol_index', 'som'], ascending=[True, True])

    # export the table to a csv file
    all_som.to_csv('../data/zaretzki/0_all/3a4_all_som.csv', index=False)
