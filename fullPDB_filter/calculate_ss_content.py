#!/usr/bin/env python

# Usage: python script.py path_to_fullPDB pdb_list(txt)
# This script calculates the helical content (in %) to determine whether a PDB can be considered mainly helical, mainly beta or coils:
# Optional outputs configured, but must be uncommented first.

# LIBRARIES IMPORTATION:
import sys
from Bio.PDB import *
parser = PDBParser(QUIET = True)
import numpy
from math import sqrt
import pandas as pd


# ARGUMENTS:
fullPDB = sys.argv[1]
pdb_list = sys.argv[2]
counter = 0
pdb_list = open(pdb_list, "r").read().splitlines()
helices = open("alphaHelices.txt", "w")
# strands = open("betaStrands.txt", "w")
# coils = open("coils.txt", "w")
results = open("ss_annotation.txt","w")
results.write(f"pdb\talpha_helical(%)\tbeta_strand(%)\tcoil(%)\n")

# LISTS:
atoms_list = list()
# CA (carbon alpha):
ca = list(); x_ca = list(); y_ca = list(); z_ca = list()
# O (carbonyl oxigen):
oxi = list(); x_oxi = list(); y_oxi = list(); z_oxi = list()
# Coordinates from centroid of CAs:
x_centroidCA = list(); y_centroidCA = list(); z_centroidCA = list()
# Coordinates from centroid of Os:
x_centroidOxi = list(); y_centroidOxi = list(); z_centroidOxi = list()


# FUNCTIONS:
# Function to calculate centroids:
#   - Centroid calculation formula: summatory of all elements / total elements (mean)
def calculateCentroids(totalAtoms, x_coord, y_coord, z_coord):
    x_list = list(); y_list = list(); z_list = list()
    for i in range(0,len(totalAtoms)-6):
        # Note: we are calculating the centroids considering coordenates for 7 residues at a time (as helices are
        # conformed by heptads).
        x = (x_coord[i] + x_coord[i+1] + x_coord[i+2] + x_coord[i+3] + x_coord[i+4] + x_coord[i+5] + x_coord[i+6]) / 7
        x_list.append(x)
        y = (y_coord[i] + y_coord[i+1] + y_coord[i+2] + y_coord[i+3] + y_coord[i+4] + y_coord[i+5] + y_coord[i+6]) / 7
        y_list.append(y)
        z = (z_coord[i] + z_coord[i+1] + z_coord[i+2] + z_coord[i+3] + z_coord[i+4] + z_coord[i+5] + z_coord[i+6]) / 7
        z_list.append(z)
    return x_list, y_list, z_list

# Function to calculate Euclidean distances: distances between two points in three-dimensional space.
#   - Euclidean distance calculation formula: sqrt((x2 – x1)^2 + (y2 – y1)^2 + … + (zn – z1)^2)
def calculateDistances(xCA, yCA, zCA, xOxi, yOxi, zOxi):
    distances_list = list()
    array_xCA = numpy.array(xCA)
    array_yCA = numpy.array(yCA)
    array_zCA = numpy.array(zCA)
    array_xOxi = numpy.array(xOxi)
    array_yOxi = numpy.array(yOxi)
    array_zOxi = numpy.array(zOxi)
    for i in range(0, len(array_xCA)):
        distance = sqrt((array_xOxi[i] - array_xCA[i])**2 + (array_yOxi[i] - array_yCA[i])**2 + (array_zOxi[i] - array_zCA[i])**2)
        distances_list.append(distance)
    return distances_list

# Function to classify structures into alpha helices, beta strands or coils and, then, calculate the helical content:
#   - Parameters for Euclidean distances:
#       - alpha helix: 2.20 Å (+-0.18 Å) ////// used: 1.72 <= ele <= 2.72 -> 2.22 sd 0.06, instead sd +-0.5
#       - beta strand: 1.39 Å (+-0.24 Å)
#       - coil: else
#   - A structure will be considered mainly helical if it contains a minimum of a 75% of helices.
def calculateHelicalContent(distances_list, pdb):
    ss_annotation = list()
    for distance in distances_list:
        if (1.7) <= distance <= (2.7):
            ss_annotation.append("H")
        elif (1.39-0.24) <= distance <= (1.39+0.24):
            ss_annotation.append("E")
        else:
            ss_annotation.append("C")
    helical_content = ss_annotation.count("H") / len(ss_annotation)
    beta_content = ss_annotation.count("E") / len(ss_annotation)
    coil_content = ss_annotation.count("C") / len(ss_annotation)
    results.write(f"{pdb}\t{round(helical_content*100,2)}\t{round(beta_content*100,2)}\t{round(coil_content*100,2)}\n")
    if helical_content >= 0.75:
        helices.write(f"{pdb}\n")
        print(f"{pdb} is considered mainly helical with a {round(helical_content*100,2)}% helical content.")
    # elif beta_content >= 0.75:
    #     strands.write(f"{pdb}\t{round(beta_content*100,2)}\n")
    #     print(f"{pdb} is considered a mainly beta with a {round(beta_content * 100, 2)}% beta content.")
    # else:
    #     coils.write(f"{pdb}\tH:{round(helical_content*100,2)}\t E:{round(beta_content*100,2)}\n")
    #     print(f"{pdb} is considered a coil with a {round(helical_content*100,2)}%/{round(beta_content * 100, 2)}% helical/beta content.")


# RUNNING SCRIPT
file = open("pre_ss_annotation.txt","r")
df = pd.read_csv(file, sep='\t', header=(0))
listrepeat = df.iloc[:,0].values.tolist()
# Iteration over a PDB list:
for pdb in pdb_list:
    if pdb not in listrepeat:
        print("NEW PDB")
        counter += 1
        ent_file = f"{fullPDB}pdb{pdb}.ent"
        print(ent_file)
        try:
            structure = parser.get_structure(pdb, ent_file)
            for model in structure:
                for chain in model:
                    for residue in chain:
                        atoms = residue.get_atoms()
                        atoms_list = {atom.get_name() for atom in atoms}
                        if all(name in atoms_list for name in ["N", "CA", "C", "O"]):
                            ca.append(residue["CA"])
                            x, y, z = residue["CA"].get_coord()
                            x_ca.append(x)
                            y_ca.append(y)
                            z_ca.append(z)
                            oxi.append(residue["O"])
                            x, y, z = residue["O"].get_coord()
                            x_oxi.append(x)
                            y_oxi.append(y)
                            z_oxi.append(z)
                        atoms_list.clear()
            x_centroidCA, y_centroidCA, z_centroidCA = calculateCentroids(ca, x_ca, y_ca, z_ca)
            x_centroidOxi, y_centroidOxi, z_centroidOxi = calculateCentroids(oxi, x_oxi, y_oxi, z_oxi)
            euclidean_distances = calculateDistances(x_centroidCA, y_centroidCA, z_centroidCA, x_centroidOxi, y_centroidOxi, z_centroidOxi)
            calculateHelicalContent(euclidean_distances, pdb)
            print(f"{pdb} has been classified. Structure number: {counter}.")
        except:
            print(f"Error: {pdb} PDB not found! Located in line {counter}.")

# Closing used files:
helices.close()
# strands.close()
# coils.close()
results.close()


