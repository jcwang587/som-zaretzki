# Description: Utility functions for the project
# The codes need to be run with schrodinger API 2024-3


import os
import shutil
import tempfile
import pandas as pd
import networkx as nx
from io import StringIO

from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem

from schrodinger.structure import StructureReader
from schrodinger.structutils.build import delete_hydrogens, add_hydrogens


def get_atom_environment(mol, atom_idx):
    """
    Get the environment of an atom (its neighbors and bond types).

    :param mol: RDKit molecule object
    :param atom_idx: Index of the atom in the molecule
    :return: Set of (neighbor_atom_symbol, bond_order) tuples
    """
    atom = mol.GetAtomWithIdx(atom_idx)

    environment_1st_shell = set()
    environment_1st_shell_list = []
    environment_2nd_shell = set()
    environment_2nd_shell_list = []
    environment_3rd_shell = set()
    environment_3rd_shell_list = []
    environment_4th_shell = set()
    environment_4th_shell_list = []
    environment_5th_shell = set()
    environment_5th_shell_list = []

    # first shell environment
    for bond_1st_shell in atom.GetBonds():
        atom_1st_shell = bond_1st_shell.GetOtherAtom(atom)
        environment_1st_shell.add(atom_1st_shell.GetSymbol())
        environment_1st_shell_list.append(atom_1st_shell.GetSymbol())
        # second shell environment
        for bond_2nd_shell in atom_1st_shell.GetBonds():
            atom_2nd_shell = bond_2nd_shell.GetOtherAtom(atom_1st_shell)
            environment_2nd_shell.add(atom_2nd_shell.GetSymbol())
            environment_2nd_shell_list.append(atom_2nd_shell.GetSymbol())
            # third shell environment
            for bond_3rd_shell in atom_2nd_shell.GetBonds():
                atom_3rd_shell = bond_3rd_shell.GetOtherAtom(atom_2nd_shell)
                environment_3rd_shell.add(atom_3rd_shell.GetSymbol())
                environment_3rd_shell_list.append(atom_3rd_shell.GetSymbol())
                # fourth shell environment
                for bond_4th_shell in atom_3rd_shell.GetBonds():
                    atom_4th_shell = bond_4th_shell.GetOtherAtom(atom_3rd_shell)
                    environment_4th_shell.add(atom_4th_shell.GetSymbol())
                    environment_4th_shell_list.append(atom_4th_shell.GetSymbol())
                    # fifth shell environment
                    for bond_5th_shell in atom_4th_shell.GetBonds():
                        atom_5th_shell = bond_5th_shell.GetOtherAtom(atom_4th_shell)
                        environment_5th_shell.add(atom_5th_shell.GetSymbol())
                        environment_5th_shell_list.append(atom_5th_shell.GetSymbol())

    # get the number of each element in the first shell
    environment_1st_shell_dict = {
        i: environment_1st_shell_list.count(i) for i in environment_1st_shell
    }
    environment_1st_shell_count = list(environment_1st_shell_dict.values())
    environment_2nd_shell_dict = {
        i: environment_2nd_shell_list.count(i) for i in environment_2nd_shell
    }
    environment_2nd_shell_count = list(environment_2nd_shell_dict.values())
    environment_3rd_shell_dict = {
        i: environment_3rd_shell_list.count(i) for i in environment_3rd_shell
    }
    environment_3rd_shell_count = list(environment_3rd_shell_dict.values())
    environment_4th_shell_dict = {
        i: environment_4th_shell_list.count(i) for i in environment_4th_shell
    }
    environment_4th_shell_count = list(environment_4th_shell_dict.values())
    environment_5th_shell_dict = {
        i: environment_5th_shell_list.count(i) for i in environment_5th_shell
    }
    environment_5th_shell_count = list(environment_5th_shell_dict.values())

    env_1 = [environment_1st_shell, environment_1st_shell_count]
    env_2 = [environment_2nd_shell, environment_2nd_shell_count]
    env_3 = [environment_3rd_shell, environment_3rd_shell_count]
    env_4 = [environment_4th_shell, environment_4th_shell_count]
    env_5 = [environment_5th_shell, environment_5th_shell_count]

    env = [env_1, env_2, env_3, env_4, env_5]

    return env


def mol2struct(rdkit_mol):
    """
    Convert RDKit molecule object to Schrodinger structure object

    :param rdkit_mol: RDKit molecule object
    :return: Schrodinger structure object
    """
    sio = StringIO()
    temp = tempfile.NamedTemporaryFile(delete=False, suffix=".sdf")

    try:
        # Export the reactant molecule to SDF format
        w = Chem.SDWriter(temp.name)
        w.write(rdkit_mol)
        w.close()

        schrodinger_structure = StructureReader.read(temp.name)
    finally:
        os.unlink(temp.name)

    return schrodinger_structure


def renew_folder(folder_path):
    """
    Remove and recreate a folder

    :param folder_path: Path to the folder
    """
    shutil.rmtree(folder_path, ignore_errors=True)
    os.makedirs(folder_path, exist_ok=True)


def metric_acc(pred, actual):
    """
    Calculate the accuracy metric.
    
    Parameters
    ----------
    pred : dict with list values
    actual : dict with list values

    Returns
    -------
    float
    """
    correct = 0
    fail_key = []
    total = len(pred)
    for idx in pred.keys():
        # check if any values in the predicted list are in the actual list
        if any(x in actual[idx] for x in pred[idx]):
            correct += 1
        # if none of the values in the predicted list are in the actual list
        else:
            fail_key.append(idx)
        
    acc = correct / total
 
    return acc, fail_key


def metric_se(pred, actual):
    """
    Calculate the sensitivity metric, which is the true positive rate.
    
    Parameters
    ----------
    pred : dict
    actual : dict
    
    Returns
    -------
    float
    """
    tp = 0
    fn = 0
    for idx in pred.keys():
        if pred[idx] in actual[idx]:
            tp += 1
        else:
            fn += 1
    se = tp / (tp + fn)
    return se


def metric_sp(pred, actual):
    """
    Calculate the specificity metric, which is the true negative rate.
    
    Parameters
    ----------
    pred : dict
    actual : dict
    
    Returns
    -------
    float
    """
    tn = 0
    fp = 0
    for idx in pred.keys():
        if pred[idx] in actual[idx]:
            tn += 1
        else:
            fp += 1
    sp = tn / (tn + fp)
    return sp


def metric_bacc(pred, actual):
    """
    Calculate the balanced accuracy metric.
    
    Parameters
    ----------
    pred : dict
    actual : dict
    
    Returns
    -------
    float
    """
    se = metric_se(pred, actual)
    sp = metric_sp(pred, actual)
    bacc = (se + sp) / 2
    return bacc


def get_3a4_reactivity(input_structure, top_k=1, export_mae=False):
    """
    Get the indices of the top k atoms with the highest reactivity.

    Parameters
    ----------
    input_structure : schrodinger.structure.StructureReader
    top_k : int, optional
    export_mae : bool, optional

    """
    delete_hydrogens(input_structure)

    property_values = pd.DataFrame(columns=['r_cyp_CYP_reactivity'])
    
    atoms = input_structure.atom
    for i, atom in enumerate(atoms):
        if 'r_cyp_CYP_reactivity' in atom.property:
            property_values.loc[i, 'r_cyp_CYP_reactivity'] = atom.property['r_cyp_CYP_reactivity'] * -1
        else:
            property_values.loc[i, 'r_cyp_CYP_reactivity'] = None

    if export_mae:
        input_structure.write(export_mae)

    property_values['r_cyp_CYP_reactivity'] = property_values['r_cyp_CYP_reactivity'].astype(float)
    top_reactivity_indices = property_values.nlargest(top_k, 'r_cyp_CYP_reactivity').index + 1

    return top_reactivity_indices.tolist()


def mol_to_nx(mol):
    """
    Convert an RDKit molecule to a NetworkX graph.

    :param mol: RDKit molecule object
    :return: NetworkX graph
    """
    g = nx.Graph()
    for atom in mol.GetAtoms():
        g.add_node(atom.GetIdx(), symbol=atom.GetSymbol(), atom=atom)
    for bond in mol.GetBonds():
        g.add_edge(
            bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), order=bond.GetBondType()
        )
    return g


def compare_molecules_entry(reactant_mol, product_mol, reaction_key):
    """
    Compare the reactant and product molecules and classify the reaction type.

    :param reactant_mol: RDKit molecule object (reactant)
    :param product_mol: RDKit molecule object (product)
    :param reaction_key: Reaction key from BioTransformer output
    """
    if reaction_key == "methylation":
        return compare_molecules_add_one_atom(reactant_mol, product_mol, "C")
    elif reaction_key == "hydroxylation":
        return compare_molecules_add_one_atom(reactant_mol, product_mol, "O")
    elif reaction_key == "oxydation":
        return compare_molecules_add_one_atom(reactant_mol, product_mol, "O")


def compare_molecules_add_one_atom(reactant_mol, product_mol, element):
    """
    Compare the reactant and product molecules and identify one new carbon atom and its bonds.

    :param reactant_mol: RDKit molecule object (reactant)
    :param product_mol: RDKit molecule object (product)
    :param element: Element to compare (e.g., 'C' for carbon)
    :return: Tuple of new atom indices and new bond tuples (atom indices) from the product structure
    """
    new_atoms = []
    new_bonds = []
    diff_score_list = []

    # Compare atoms between reactant and product by environment (atom type and bond connections)
    for atom_idx in range(product_mol.GetNumAtoms()):
        diff_score = 0

        product_atom = product_mol.GetAtomWithIdx(atom_idx)
        if product_atom.GetSymbol() == element:
            # Check if this carbon exists in the reactant by comparing environments
            product_env = get_atom_environment(product_mol, atom_idx)[:5]
            product_env_1, product_env_2, product_env_3, product_env_4, product_env_5 = product_env

            for reactant_idx in range(reactant_mol.GetNumAtoms()):
                if (
                    product_atom.GetSymbol()
                    == reactant_mol.GetAtomWithIdx(reactant_idx).GetSymbol()
                ):
                    reactant_env = (get_atom_environment(reactant_mol, reactant_idx))[:5]
                    reactant_env_1, reactant_env_2, reactant_env_3, reactant_env_4, reactant_env_5 = reactant_env
                    # Compare the first and second shell environments
                    if product_env_1 != reactant_env_1:
                        diff_score += 10000
                    if product_env_2 != reactant_env_2:
                        diff_score += 1000
                    if product_env_3 != reactant_env_3:
                        diff_score += 100
                    if product_env_4 != reactant_env_4:
                        diff_score += 10
                    if product_env_5 != reactant_env_5:
                        diff_score += 1
            diff_score_list.append(diff_score)
            new_atoms.append(atom_idx)

    # get the new atoms which has the largest difference score
    new_atoms_idx = (diff_score_list.index(max(diff_score_list)))
    new_atoms = [new_atoms[new_atoms_idx]]

    # Find new bonds involving the new carbon atom
    for atom_idx in new_atoms:
        for bond in product_mol.GetAtomWithIdx(atom_idx).GetBonds():
            begin_idx = bond.GetBeginAtomIdx()
            end_idx = bond.GetEndAtomIdx()
            bond_order = bond.GetBondType()
            if begin_idx in new_atoms or end_idx in new_atoms:
                new_bonds.append((begin_idx, end_idx, {"order": bond_order}))

    return new_atoms, new_bonds


def highlight_save_structure_from_product(
    product_mol,
    new_atoms,
    new_bonds,
    export_driver,
    reactant_mol=None,
    draw_svg=True,
    output_filename="highlighted_structure",
):
    """
    Highlight the new atoms and bonds in the molecule and save the structure.
    Two modes are availabe, (1) one is the single job mode which can directly 
    export the sdf file through rdkit or mae file of the product structure through schrodinger,
    (2) another one is the batch mode which can return the reaction site index for the reactant molecule.

    :param product_mol: RDKit molecule object (product)
    :param new_atoms: Set of new atom indices
    :param new_bonds: Set of new bond tuples (atom indices)
    :param reactant_mol: RDKit molecule object (reactant)
    :param export_driver: Drawing driver ('rdkit' or 'schrodinger' or 'none')
    :param draw_svg: Export the structure as an SVG image (True or False)
    :param output_filename: File name to save the output (e.g., 'output.sdf', 'output.png')
    """

    # Convert the new edges to bond indices (RDKit uses bond indices, not atom indices for highlighting)
    bond_indices = []
    for bond in product_mol.GetBonds():
        if (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()) in [
            (b[0], b[1]) for b in new_bonds
        ] or (bond.GetEndAtomIdx(), bond.GetBeginAtomIdx()) in [
            (b[0], b[1]) for b in new_bonds
        ]:
            bond_indices.append(bond.GetIdx())

    # Highlight atoms and bonds in the structure (new_atoms for atoms and bond_indices for bonds)
    bond_atoms = list(
        set([bond[0] for bond in new_bonds] + [bond[1] for bond in new_bonds])
    )
    old_atoms = list(set(bond_atoms) - set(new_atoms))

    # Save the SVG image
    if draw_svg:
        highlight = list(set(new_atoms + old_atoms))
        highlight_colors = {
            atom: (0, 1, 0) for atom in old_atoms
        }  # Green for old atoms, red for new atoms

        drawer = Draw.rdMolDraw2D.MolDraw2DSVG(500, 300)
        drawer.DrawMolecule(
            product_mol,
            highlightAtoms=highlight,
            highlightBonds=bond_indices,
            highlightAtomColors=highlight_colors,
        )
        drawer.FinishDrawing()

        svg = drawer.GetDrawingText()
        with open("product_" + output_filename + ".svg", "w") as f:
            f.write(svg)

        print(f"Structure saved as {output_filename}.svg")

    if export_driver == "none":
        # Remove the new atoms from the product molecule
        product_mol = Chem.RWMol(product_mol)
        for atom_idx in new_atoms:
            product_mol.RemoveAtom(atom_idx)

        new_atoms_product = new_atoms.copy()
        old_atoms_product = old_atoms.copy()

        # Find the exact same atoms in the reactant as the old atoms in the product
        if new_atoms_product[0] < old_atoms_product[0]:
            print("Warning: New atom index is smaller than the old atom index")
            old_atoms_product[0] -= 1

        print(f"Old atom index: {old_atoms_product[0]}")

        old_atom_product_env = get_atom_environment(product_mol, old_atoms_product[0])[:5]

        old_atoms_reactant = []
        for reactant_idx in range(reactant_mol.GetNumAtoms()):
            reactant_env = get_atom_environment(reactant_mol, reactant_idx)[:5]
            if old_atom_product_env == reactant_env:
                old_atoms_reactant.append(reactant_idx)

        # Check if old_atoms_reactant is empty, if so error message
        if not old_atoms_reactant:
            print("Old atom not found in the reactant molecule environment")

        return old_atoms_reactant

    elif export_driver == "rdkit":
        reactive_list = [
            1 if atom.GetIdx() in old_atoms else 0 for atom in product_mol.GetAtoms()
        ]

        # Add reactive property to the new atoms
        product_mol = Chem.RWMol(product_mol)
        for atom_idx in range(len(reactive_list)):
            product_mol.GetAtomWithIdx(atom_idx).SetIntProp(
                "reactive_site", reactive_list[atom_idx]
            )

        Chem.CreateAtomIntPropertyList(product_mol, "reactive_site")

        # Remove the new atoms from the molecule
        for atom_idx in new_atoms:
            product_mol.RemoveAtom(atom_idx)

        AllChem.EmbedMolecule(product_mol)

        w = Chem.SDWriter(output_filename + ".sdf")
        w.write(product_mol)
        w.close()

        print(f"Structure saved as {output_filename}.sdf")

    elif export_driver == "schrodinger":
        # Convert the RDKit molecule to a Schrodinger structure
        AllChem.EmbedMolecule(product_mol)

        schrodinger_structure = mol2struct(product_mol)
        reactive_bool_list = [
            True if atom.GetIdx() in old_atoms else False
            for atom in product_mol.GetAtoms()
        ]

        # Add reactive property to the atoms
        for atom_idx in range(len(reactive_bool_list)):
            schrodinger_structure.atom[atom_idx + 1].property[
                "b_user_oxidation_biotransformer"
            ] = reactive_bool_list[atom_idx]

        # Remove the new atoms from the molecule
        new_atoms_slice = slice(min(new_atoms), max(new_atoms) + 1)
        new_atoms_schrondinger = list(schrodinger_structure.atom)[new_atoms_slice]
        for atom_idx in new_atoms:
            schrodinger_structure.deleteAtoms(new_atoms_schrondinger)

        # Add hydorgen atoms to the structure
        add_hydrogens(schrodinger_structure)

        # Export the structure with the reactive property to Maestro format
        schrodinger_structure.write(output_filename + ".mae")


def biotransformer_oxidation_site(input_csv, output_file):
    """
    Load the BioTransformer output CSV file and identify the oxidation sites in the product molecules.

    :param input_csv: Path to the BioTransformer output CSV file
    :param output_file: Path to save the output Maestro file for property assignment or mol2 file for gold docking
    """
    # Load the BioTransformer output CSV file
    biotransformer_output = pd.read_csv(input_csv)

    # Filter the reactions by the keyword CYP enzyme and oxidation reaction
    enzymes = biotransformer_output["Enzyme(s)"]
    cyp_enzymes = [i for i in range(len(enzymes)) if "CYP" in enzymes[i]]
    biotransformer_output = biotransformer_output.iloc[cyp_enzymes]
    biotransformer_output = biotransformer_output.reset_index(drop=True)
    reaction = biotransformer_output["Reaction"]
    oxidation_key_reactions = [
        i for i in range(len(reaction)) if "oxidation" in reaction[i].lower()
    ]
    hydroxyl_key_reactions = [
        i for i in range(len(reaction)) if "hydroxylation" in reaction[i].lower()
    ]
    epoxidation_key_reactions = [
        i for i in range(len(reaction)) if "epoxidation" in reaction[i].lower()
    ]
    oxidation_reactions = list(
        set(oxidation_key_reactions + hydroxyl_key_reactions)
        - set(epoxidation_key_reactions)
    )

    # check if the oxidation_reactions is empty, if so error message and return
    if not oxidation_reactions:
        print("No oxidation reactions found in the BioTransformer output")
        return

    biotransformer_output = biotransformer_output.iloc[oxidation_reactions]
    biotransformer_output = biotransformer_output.reset_index(drop=True)

    reactant_smiles = biotransformer_output["Precursor SMILES"]
    product_smiles = biotransformer_output["SMILES"]

    oxidation_sites = []

    # Convert SMILES to RDKit Mol objects (without adding explicit hydrogens)
    for i in range(len(reactant_smiles)):
        reactant_mol = Chem.MolFromSmiles(reactant_smiles[i])
        product_mol = Chem.MolFromSmiles(product_smiles[i])

        # Compare the molecules specifically for the reaction type
        new_atoms, new_bonds = compare_molecules_entry(
            reactant_mol, product_mol, "oxydation"
        )

        # Print the results
        print(f"New atoms (indices) in the product: {new_atoms}")
        print(f"New bonds in the product: {new_bonds}")

        # Highlight and save the structure with the new methyl group
        oxydation_site = highlight_save_structure_from_product(
            product_mol=product_mol,
            reactant_mol=reactant_mol,
            new_atoms=new_atoms,
            new_bonds=new_bonds,
            export_driver="none",
            draw_svg=False,
        )
        oxidation_sites.append(oxydation_site[0])

    # Print the count of the True values in the oxidation
    print(oxidation_sites)
    print(set(oxidation_sites))

    AllChem.EmbedMolecule(reactant_mol)

    # convert the rdkit molecule to schrodinger structure
    schrodinger_structure = mol2struct(reactant_mol)

    reactive_bool_list = [
        1 if atom.GetIdx() in oxidation_sites else 0
        for atom in reactant_mol.GetAtoms()
    ]

    # Add reactive property to the atoms
    print(len(reactive_bool_list))
    print(reactive_bool_list)
    for atom_idx in range(len(reactive_bool_list)):
        schrodinger_structure.atom[atom_idx + 1].property[
            "r_user_oxidation_biotransformer"
        ] = reactive_bool_list[atom_idx]

    # Export the structure with the reactive property to Maestro format or mol2
    schrodinger_structure.write(output_file)


def get_primary_som(structure):
    """
    Get the primary SOM of the molecule

    :param structure: Schrodinger structure object
    :return: List of primary SOM indices
    """
    # Get the SOM of the molecule
    primary_som_all = []
    # Iterate through the structures and get the PRIMARY_SOM property
    primary_som_int = structure.property.get('i_sd_PRIMARY\_SOM', 'N/A')
    if primary_som_int == 'N/A':
        primary_som_string = structure.property.get('s_sd_PRIMARY\_SOM', 'N/A')
        primary_som_list = list(map(int, primary_som_string.split()))
        primary_som_all.append(primary_som_list)
    else:       
        primary_som_all.append([primary_som_int])

    return primary_som_all[0]

def get_secondary_som(structure):
    """
    Get the secondary SOM of the molecule

    :param structure: Schrodinger structure object
    :return: List of secondary SOM indices
    """
    # Get the SOM of the molecule
    secondary_som_all = []

    # Check if the SECONDARY_SOM property exists
    if 'i_sd_SECONDARY\_SOM' in structure.property.keys() or 's_sd_SECONDARY\_SOM' in structure.property.keys():
        secondary_som_int = structure.property.get('i_sd_SECONDARY\_SOM', 'N/A')
        if secondary_som_int == 'N/A':
            secondary_som_string = structure.property.get('s_sd_SECONDARY\_SOM', 'N/A')
            secondary_som_list = list(map(int, secondary_som_string.split()))
            secondary_som_all.append(secondary_som_list)
        else:
            secondary_som_all.append([secondary_som_int])
        return secondary_som_all[0]
    else:
        return []


def get_tertiary_som(structure):
    """
    Get the tertiary SOM of the molecule

    :param structure: Schrodinger structure object
    :return: List of tertiary SOM indices
    """
    # Get the SOM of the molecule
    tertiary_som_all = []

    # Check if the TERTIARY_SOM property exists
    if 'i_sd_TERTIARY\_SOM' in structure.property.keys() or 's_sd_TERTIARY\_SOM' in structure.property.keys():
        tertiary_som_int = structure.property.get('i_sd_TERTIARY\_SOM', 'N/A')
        if tertiary_som_int == 'N/A':
            tertiary_som_string = structure.property.get('s_sd_TERTIARY\_SOM', 'N/A')
            tertiary_som_list = list(map(int, tertiary_som_string.split()))
            tertiary_som_all.append(tertiary_som_list)
        else:
            tertiary_som_all.append([tertiary_som_int])
        return tertiary_som_all[0]  
    else:
        return []

