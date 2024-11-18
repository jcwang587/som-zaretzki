from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
import cairosvg
import shutil

import sys

sys.path.append("../src")
from utils import renew_folder

renew_folder("../data/zaretzki/0_all/svg")
renew_folder("../data/zaretzki/0_all/png")

# Load sdf file with rdkit
sdf_file = "../data/zaretzki/0_all/3A4.sdf"
suppl = Chem.SDMolSupplier(sdf_file)

# Iterate over each molecule in the SDF file
for index, mol in enumerate(suppl):
    if mol is None:
        continue  # Skip if molecule is None

    # Compute 2D coordinates
    AllChem.Compute2DCoords(mol)

    # Initialize lists for SOM atom indices
    primary_som_atoms = []
    secondary_som_atoms = []
    tertiary_som_atoms = []

    # Get PRIMARY_SOM atom indices
    if mol.HasProp("PRIMARY_SOM"):
        primary_som = mol.GetProp("PRIMARY_SOM")
        primary_som_atoms = [int(i) - 1 for i in primary_som.strip().split()]

    # Get SECONDARY_SOM atom indices
    if mol.HasProp("SECONDARY_SOM"):
        secondary_som = mol.GetProp("SECONDARY_SOM")
        secondary_som_atoms = [int(i) - 1 for i in secondary_som.strip().split()]

    # Get TERTIARY_SOM atom indices
    if mol.HasProp("TERTIARY_SOM"):
        tertiary_som = mol.GetProp("TERTIARY_SOM")
        tertiary_som_atoms = [int(i) - 1 for i in tertiary_som.strip().split()]

    # Combine all SOM atom indices
    highlight_atoms = primary_som_atoms + secondary_som_atoms + tertiary_som_atoms

    # Define colors for each SOM type
    primary_color = (0.7, 0, 0)
    secondary_color = (0, 0.7, 0)
    tertiary_color = (0, 0.7, 0.7)  # cyan

    # Create a color dictionary mapping atom indices to colors
    highlight_colors = {}
    for idx in primary_som_atoms:
        highlight_colors[idx] = primary_color
    for idx in secondary_som_atoms:
        highlight_colors[idx] = secondary_color
    for idx in tertiary_som_atoms:
        highlight_colors[idx] = tertiary_color

    # Draw the molecule with highlighted atoms without bonds
    drawer = Draw.rdMolDraw2D.MolDraw2DSVG(500, 300)
    opts = drawer.drawOptions()
    for i in range(mol.GetNumAtoms()):
        opts.atomLabels[i] = f"{mol.GetAtomWithIdx(i).GetSymbol()}{i+1}"
    drawer.DrawMolecule(
        mol,
        highlightAtoms=highlight_atoms,
        highlightAtomColors=highlight_colors,
        highlightBonds=[],
    )
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()

    # Save the SVG image
    mol_name = mol.GetProp("_Name") if mol.HasProp("_Name") else "molecule"
    with open(
        f"../data/zaretzki/0_all/svg/{index+1:03d}_{mol_name}_som.svg", "w"
    ) as f:
        f.write(svg)

    cairosvg.svg2png(
        url=f"../data/zaretzki/0_all/svg/{index+1:03d}_{mol_name}_som.svg",
        write_to=f"../data/zaretzki/0_all/png/{index+1:03d}_{mol_name}_som.png",
    )


# remove the svg folder
shutil.rmtree("../data/zaretzki/0_all/svg")
