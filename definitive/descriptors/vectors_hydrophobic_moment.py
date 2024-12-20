# !usr/bin/env python

# Usage: python script.py fullPDB_path PDB_list(.txt)

# LIBRARIES IMPORTATION:
from Bio.PDB import *
parser = PDBParser(QUIET=True)
io = PDBIO()
from Bio.SeqUtils import IsoelectricPoint
sr = ShrakeRupley()
import traceback
import os, sys, subprocess
import pandas as pd
import math
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

# SCALES:
consensus1982 = {
    "ALA": 0.25,
    "CYS": 0.04,
    "ASP": -0.72,
    "GLU": -0.62,
    "PHE": 0.61,
    "GLY": 0.16,
    "HIS": -0.40,
    "ILE": 0.73,
    "LYS": -1.10,
    "LEU": 0.53,
    "MET": 0.26,
    "ASN": -0.64,
    "PRO": -0.07,
    "GLN": -0.69,
    "ARG": -1.76, # in the original paper the value is rounded to -1.8
    "SER": -0.26,
    "THR": -0.18,
    "VAL": 0.54,
    "TRP": 0.37,
    "TYR": 0.02
}

hydrophobic = ["ALA", "VAL", "LEU", "ILE", "PHE", "TRP", "MET", "PRO","CYS", "GLY", "MSE"]

# ARGUMENTS:
fullPDB = sys.argv[1]
df = sys.argv[2]
df = pd.read_csv(df, sep="\t")
dir = os.getcwd()
errors = open("errors.csv","w")
if "graphical_output" not in os.listdir(dir):
    subprocess.run(f"mkdir graphical_output", shell=True)
if "single_chain_PDBs" not in os.listdir(dir):
    subprocess.run(f"mkdir single_chain_PDBs", shell=True)

# FUNCTIONS:
# Function to calculate Euclidean distances: distances between two points in three-dimensional space.
#   - Euclidean distance calculation formula: sqrt((x2 – x1)^2 + (y2 – y1)^2 + … + (zn – z1)^2)
def calculateDistances(xCA, yCA, zCA, xCB, yCB, zCB, xOxi, yOxi, zOxi, residues, pdb,oligomer):
    distances_list21 = list()
    distances_list31 = list()
    # Prepare lists for quiver plot (to accumulate all displacement vectors)
    vectors_start = []
    vectors_start_centroids = []
    vectors_displacement21 = []
    vectors_displacement31 = []
    colors = []
    # Loop over the coordinates to calculate distances and visualize vectors
    for i in range(len(xCA)-1):
        # Define the vectors for the points at index 'i'
        vector1 = np.array([xCA[i], yCA[i], zCA[i]])
        vector2 = np.array([xCB[i], yCB[i], zCB[i]])
        #vector3 = np.array([xOxi[i], yOxi[i], zOxi[i]])
        vector3 = np.array([xCA[i+1], yCA[i+1], zCA[i+1]])
        # Calculate the Euclidean distance between vector1 and vector2
        distance21 = np.linalg.norm(vector2 - vector1)
        distances_list21.append(distance21)
        # Calculate the Euclidean distance between vector1 and vector3
        distance31 = np.linalg.norm(vector3 - vector1)
        distances_list31.append(distance31)
        # Create the displacement vector (difference between vector2 and vector1)
        displacement21 = vector2 - vector1
        # Create the displacement vector (difference between vector3 and vector1)
        displacement31 = vector3 - vector1
        # Accumulate vectors for quiver plot
        vectors_start.append(vector1)
        vectors_displacement21.append(displacement21)
        vectors_displacement31.append(displacement31)
        if residues[i] in hydrophobic:
            colors.append("b")
        else:
            colors.append("r")
    # Create a 3D plot to visualize the vectors
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    # Plot all displacement vectors in a single quiver plot
    num = 0
    for start, displacement in zip(vectors_start, vectors_displacement21):
        ax.quiver(start[0], start[1], start[2], displacement[0], displacement[1], displacement[2],
                  pivot='tail', length=1, color=colors[num], arrow_length_ratio=0.5)
        num += 1
    for start, displacement in zip(vectors_start, vectors_displacement31):
        ax.quiver(start[0], start[1], start[2], displacement[0], displacement[1], displacement[2],
                  pivot='tail', length=1, color="g", arrow_length_ratio=0.1)
    # Set the limits for the axes based on your data range
    ax.set_xlim([min(xCA + xCB) - 2, max(xCA + xCB) + 2])
    ax.set_ylim([min(yCA + yCB) - 2, max(yCA + yCB) + 2])
    ax.set_zlim([min(zCA + zCB) - 2, max(zCA + zCB) + 2])
    # Create custom legend handles
    legend = [
        Line2D([0], [0], color="b", lw=4, label="Hydrophobic"),
        Line2D([0], [0], color="r", lw=4, label="Other"),
        Line2D([0], [0], color="g", lw=4, label="Backbone")
    ]
    # Show the plot with the custom legend
    ax.legend(handles=legend, loc="best", title=f"{pdb.upper()}: {oligomer}")
    ax.set_aspect('equal', 'box')
    #plt.tight_layout()
    #plt.title(f"{pdb}", loc="center")
    plt.savefig(f"graphical_output/hydrophobicity_vectors_{pdb}.png")
    plt.show()
    return distances_list21, distances_list31

def calculateNetCharge(structure):
    IP = IsoelectricPoint.IsoelectricPoint(structure)
    charge = round(IP.charge_at_pH(7.0), 2)
    return charge


# EXECTION:
# Iteration through all the PDBs list:
counter = 0
errors_num = 0
for i in df.index:
    pdb = df.loc[i,"pdb"]
    oligomer = df.loc[i,"oligomer"]
    counter += 1
    residues_list = []
    try:
        ent_file = os.path.join(fullPDB,f"pdb{pdb}.ent")
        structure = parser.get_structure(f"{pdb}", f"{ent_file}")
        try:
            for model in structure:
                for chain in model:
                    if Chain.Chain.get_id(chain) == "A":
                        # CALCULATE NET CHARGE FROM SINGLE CHAIN:
                        io.set_structure(chain)
                        io.save(f"single_chain_PDBs/pdb{pdb}_A.pdb")
                        df.at[i, "net_charge"] = calculateNetCharge(f"{dir}/single_chain_PDBs/pdb{pdb}_A.pdb")
        except:
            errors_num += 1
            errors.write(f"{pdb}\n")
            traceback.print_exc()
            print(f"Error in {pdb.upper()}. Could not process structure number {counter}.")
        try:
            for model in structure:
                for chain in model:
                    if Chain.Chain.get_id(chain) == "A":
                        # LISTS:
                        x_CA = []
                        y_CA = []
                        z_CA = []
                        x_CB = []
                        y_CB = []
                        z_CB = []
                        x_Oxi = []
                        y_Oxi = []
                        z_Oxi = []
                        for residue in chain:
                            atoms = residue.get_atoms()
                            atoms_list = {atom.get_name() for atom in atoms}
                            if all(name in atoms_list for name in ["N", "CA", "C", "O"]):
                                residues_list.append(residue.get_resname())
                                x, y, z = residue["CA"].get_coord()
                                x_CA.append(x)
                                y_CA.append(y)
                                z_CA.append(z)
                                if residue.get_resname() != "GLY":
                                    x, y, z = residue["CB"].get_coord()
                                    x_CB.append(x)
                                    y_CB.append(y)
                                    z_CB.append(z)
                                else:
                                    x_CB.append(x)
                                    y_CB.append(y)
                                    z_CB.append(z)
                                x, y, z = residue["O"].get_coord()
                                x_Oxi.append(x)
                                y_Oxi.append(y)
                                z_Oxi.append(z)
                        distances21, distances31 = calculateDistances(x_CA,y_CA,z_CA,x_CB,y_CB,z_CB,x_Oxi,y_Oxi,z_Oxi,residues_list,pdb,oligomer)
                        atoms_list.clear()
            print(f"{pdb.upper()} has been correctly processed. Structure number {counter}.")
        except:
            errors_num += 1
            errors.write(f"{pdb}\n")
            traceback.print_exc()
            print(f"Error in {pdb.upper()}. Could not process structure number {counter}.")
    except:
        errors_num += 1
        errors.write(f"{pdb}\n")
        traceback.print_exc()
        print(f"Error in {pdb.upper()}. Could not process structure number {counter}.")

print(f"\nOUTPUT SUMMARY:")
print(f"Structures correctly processed: {counter - errors_num}.")
print(f"Structures not correctly classified: {errors_num}.")
errors.close()
df.to_csv(f"output_descriptors.csv", index=False, sep="\t")



