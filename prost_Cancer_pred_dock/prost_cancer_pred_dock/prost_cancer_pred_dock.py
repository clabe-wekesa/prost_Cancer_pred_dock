#!/usr/bin/env python
# coding: utf-8

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
from padelpy import padeldescriptor
import joblib
import csv
from rdkit.Chem import AllChem
from meeko import MoleculePreparation, PDBQTWriterLegacy
from pymol import cmd
from vina import Vina
import sys
import subprocess

sys.path.insert(1, "/Jupyter_Dock/utilities")
from Jupyter_Dock.utilities.utils import getbox, pdbqt_to_sdf

# constants
docking_protein = "5gs4.pdbqt"

def calculate_lipinski_descriptors(smiles):
    """Lipinski descriptors: A set of molecular properties used to assess the drug-likeness or pharmacokinetic profile of a chemical compound

    Params
    ------
    smiles: string: An rdkit valid canonical SMILES or chemical structure a compound.

    Usage
    -----
    from prot_cancer_pred_dock import calculate_lipinski_descriptors

    calculate_lipinski_descriptors("Oc1ccc2c(c1)S[C@H](c1ccco1)[C@H](c1ccc(OCCN3CCCCC3)cc1)O2")
    """

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("You entered an invalid SMILES string")

    else:
        descriptors = {
            "Molecular Weight": Descriptors.MolWt(mol),
            "LogP": Descriptors.MolLogP(mol),
            "Num H Donors": Descriptors.NumHDonors(mol),
            "Num H Acceptors": Descriptors.NumHAcceptors(mol),
            "Num Rotatable Bonds": Descriptors.NumRotatableBonds(mol),
            "Carbon Count": Descriptors.HeavyAtomCount(mol),
            "Oxygen Count": sum(
                1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8
            ),
        }

        aliases = {
            "Molecular Weight": "Molecular Weight",
            "LogP": "LogP",
            "Num H Donors": "Number Hydrogen Bond Donors",
            "Num H Acceptors": "Number of Hydrogen Bond Acceptors",
            "Num Rotatable Bonds": "Number of Rotatable Bonds",
            "Carbon Count": "Carbon Count",
            "Oxygen Count": "Oxygen Count",
        }

        df = pd.DataFrame(
            {
                "Descriptor": list(descriptors.keys()),
                "Value": list(descriptors.values()),
            }
        )
        df["Descriptor"] = df["Descriptor"].map(aliases)

        return df


def predict_pIC50(smiles):
    """Prediction model is based on RandomForest regression constructed using a collection of all known cannonical SMILES that interact with Oestrogen Receptor alpha protein stored in the ChEMBL database.

    Params
    ------
    smiles: string: An rdkit valid canonical SMILES or chemical structure a compound.

    Usage
    -----
    from prost_cancer_pred_dock import predict_pIC50

    predict_pIC50("Oc1ccc2c(c1)S[C@H](c1ccco1)[C@H](c1ccc(OCCN3CCCCC3)cc1)O2")
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("You entered an invalid SMILES string")

    else:
        # Create a CSV file with the SMILES string and a placeholder compound name
        filename = "molecule.smi"
        data = [[smiles + "\t" + "Compound_name"]]
        with open(filename, "w", newline="") as file:
            writer = csv.writer(file)
            writer.writerows(data)

        # Process the fingerprints
        xml_files = glob.glob("fingerprints_xml/*.xml")
        xml_files.sort()
        FP_list = [
            "AtomPairs2DCount",
            "AtomPairs2D",
            "EState",
            "CDKextended",
            "CDK",
            "CDKgraphonly",
            "KlekotaRothCount",
            "KlekotaRoth",
            "MACCS",
            "PubChem",
            "SubstructureCount",
            "Substructure",
        ]
        fp = dict(zip(FP_list, xml_files))
        fingerprint = "Substructure"
        fingerprint_output_file = "".join([fingerprint, ".csv"])
        fingerprint_descriptortypes = fp[fingerprint]

        padeldescriptor(
            mol_dir="molecule.smi",
            d_file=fingerprint_output_file,
            descriptortypes=fingerprint_descriptortypes,
            detectaromaticity=True,
            standardizenitro=True,
            standardizetautomers=True,
            removesalt=True,
            log=True,
            fingerprints=True,
        )

        # Load the processed data and predict
        data = pd.read_csv("Substructure.csv")
        X = data.drop(columns=["Name"])

        loaded_model = joblib.load("padel_model.joblib")
        y_pred = loaded_model.predict(X)
        predicted_value = y_pred[0]
        predicted_value = format(predicted_value, ".2f")
        return predicted_value


def prot_lig_docking(smiles):
    """Docking procedure is performed by Autodock Vina on the Oestrogen Receptor alpha protein, pdb_id: 5gs4 prepared for docking by the Angel Moreno's Jupyter_Dock scripts.

    Params
    ------
    smiles: string, an rdkit valid canonical SMILES or chemical structure a compound.
    color: string, any matplotlib colors; default='spectrum'
    tubes: Boolean, protein visualization as cylindrical tubes, default=False

    Usage
    -----
    from prost_cancer_pred_dock import prot_lig_docking

    prot_lig_docking("Oc1ccc2c(c1)S[C@H](c1ccco1)[C@H](c1ccc(OCCN3CCCCC3)cc1)O2")
    """

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("You entered an invalid SMILES string")
    # Convert SMILES to molecule object
    else:
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=42)

        # Write the molecule to an SDF file
        writer = Chem.SDWriter("ligand_clean.sdf")
        writer.write(mol)
        writer.close()

        # Prepare the ligand
        mol_supplier = Chem.SDMolSupplier("ligand_clean.sdf", removeHs=False)
        preparator = MoleculePreparation()

        for mol in mol_supplier:
            mol_setups = preparator.prepare(mol)
            for setup in mol_setups:
                pdbqt_tuple = PDBQTWriterLegacy.write_string(setup)
                pdbqt_string = pdbqt_tuple[0]

                # Save pdbqt_string to the ligand.pdbqt file
                with open("ligand.pdbqt", "w") as pdbqt_file:
                    pdbqt_file.write(pdbqt_string)

        cmd.load(filename=docking_protein, format="pdb", object="prot")
        cmd.load(filename="ligand_clean.sdf", format="sdf", object="lig")
        center, size = getbox(selection="lig", extending=5.0, software="vina")
        cmd.delete("all")

        v = Vina(sf_name="vina")
        v.set_receptor(vina_input)
        v.set_ligand_from_file("ligand.pdbqt")
        v.compute_vina_maps(
            center=[center["center_x"], center["center_y"], center["center_z"]],
            box_size=[size["size_x"], size["size_y"], size["size_z"]]
        )

        v.dock(exhaustiveness=10, n_poses=10)
        v.write_poses("5gs4_ligand_vina_out.pdbqt", n_poses=10, overwrite=True)
        pdbqt_to_sdf(pdbqt_file="5gs4_ligand_vina_out.pdbqt", output="5gs4_ligand_vina_out.sdf")

        results = Chem.SDMolSupplier("5gs4_ligand_vina_out.sdf")

        p = Chem.MolToMolBlock(results[0], False)

        if results:
            return "Best docking score: {}".format(results[0].GetProp("Score"))
        else:
            return "No results available"

# Run PyMOL with the specified PDBQT files
def vizualize_dock_results():
    """Visualization and protein-ligand interaction in pymol. Users should only run this function after running prot_lig_docking function.

    Usage
    -----
    from prost_cancer_pred_dock import vizualize_dock_results

    vizualize_dock_results()
    """
    return subprocess.run(["pymol", docking_protein, "5gs4_ligand_vina_out.pdbqt"])       