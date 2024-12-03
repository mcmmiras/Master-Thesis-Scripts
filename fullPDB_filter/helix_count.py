#!/user/bin/env python

# Usage: python script.py fullPDB_path PDB_list option[filter/remove]

import os
import sys
import re
import pandas as pd
import subprocess
import traceback
from Bio.PDB import *
parser = PDBParser(QUIET = True)
import numpy
from math import sqrt


# ARGUMENTS:
pdb_list = sys.argv[2]
pdb_list = pd.read_csv(pdb_list, sep="\t", header=0)
fullPDB = sys.argv[1]
option = sys.argv[3]
dir = os.getcwd()
if f"{option}_helix_extractions" not in os.listdir(dir):
    subprocess.run(f"mkdir {option}_helix_extractions", shell=True)
    subprocess.run(f"mkdir {option}_helix_extractions/helices_{option}_output", shell=True)
helices_annotation = pd.DataFrame(columns=["pdb","helix_number","helix_length"])

# LISTS:
atoms_list = list()
residuesName_list = list()
residuesPosition_list = list()
resHelices_list = list()
chain_list = list()
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
def calculateHelicalContent(distances_list, resnames, respos, chains,pdb):
    counter = 0
    prev_chain = chains[0]
    helix_total = 0
    lengths = list()
    for i, pos in enumerate(respos):
        print(pos)
        try:
            if (1.7) <= distances_list[i] <= (2.7):
                print(pos)
                last_helix_res = pos
                ss = "H"
                counter += 1
                if chains[i] != prev_chain or pos != last_helix_res:
                    print("Yes")
                    counter = 1
                    prev_chain = chain[i]
                if counter == 1:
                    helix_total += 1
                    lengths.append(int(last_helix_res))
                helices.write(f"{respos[i]}\t{resnames[i]}\t{chain}\t{ss}\t{counter}\n")
            elif (1.39-0.24) <= distances_list[i] <= (1.39+0.24):
                ss = "E"
                counter = 0
                #helices.write(f"{respos[i]}\t{resnames[i]}\t{chain[i]}\t{ss}\n")

            else:
                ss = "C"
                counter = 0
                #helices.write(f"{respos[i]}\t{resnames[i]}\t{chain[i]}\t{ss}\n")

        except:
            helices.write(f"{respos[i]}\t{resnames[i]}\t{chain}\t \t \n")
    helices.write(f"\nTotal number of helices: {helix_total}.\n")
    print(lengths)
    return helix_total, lengths


# EXECUTION:
counter = 0
for i in pdb_list.index:
    pdb = pdb_list["pdb"][i]
    ent_file = os.path.join(fullPDB,f"pdb{pdb}.ent")
    helices = open(f"{option}_helix_extractions/helices_{option}_output/helices_{pdb}.txt", "w")
    counter += 1
    if option == "pdb":
        with open(ent_file, "r") as file:
            total_helices = 0
            file = file.read()
            pattern = r"HELIX..........................................................................."
            for match in re.finditer(pattern, file):
                total_helices += 1
                helices_annotation.at[total_helices-1, "helix_number"] = total_helices
                line = match.group()
                helices.write(f"{line}\n")
                pattern = r"                                ...."
                match = re.search(pattern,line)
                helix_length = match.group()
                helices_annotation.at[total_helices-1,"helix_length"] = helix_length
                helices_annotation.at[total_helices-1,"pdb"] = pdb
            helices_annotation.to_csv(f"{option}_{pdb}_helices_annotation.csv", sep="\t")
            helices_annotation.drop(helices_annotation.index, inplace=True)
            # Annotating each helix number:
            # Annotating each helix length:
            # pattern = r"... . ....  ... . .... ..                                  ..    "
            # match2 = re.search(pattern,line)
            # match2 = match2.group()
            # resS = (re.search(r"...", match2)).group()
            # resE = ((re.search(r"  \S\S\S", match2)).group()).split(" ")[1]
            # chainS = ((re.search(r" .", match2)).group()).split(" ")[1]
            # chainE = (((re.search(r" \S\S\S .", match2)).group()).split(" ")[1]).split(" ")[0]
            # resposS = (((re.search(r"  ...  ", match2)).group()).split(" ")[1]).split(" ")[0]
            # resposE = (((re.search(r"  ...  . ", match2)).group()))
            # length = ((re.search(r"                                ....", match2)).group()).split(" ")[1]
            # print(resposE)
            # helices.write(f"{resposS}\t{resS}\t{chainS}\t{resposE}\t{resE}\t{chainE}\t{length}\n")

            #helices_annotation.at[n+total_helices-1, "total_helices"] = total_helices



            # Annotating total number of helices:
            print(total_helices)
            helices.close()

            print(f"{pdb} has been annotated. Structure number: {counter}.")

    if option == "dssp":
        helices_list = 1
        matches_list = list()
        if not f"dssp_output" in os.listdir(os.path.join(dir,"dssp_helix_extractions")):
            subprocess.run(f"mkdir dssp_helix_extractions/dssp_output", shell=True)
        dssp_file = f"dssp_helix_extractions/dssp_output/{pdb}.dssp"
        if f"{pdb}.dssp" not in os.listdir(os.path.join(dir, "dssp_helix_extractions/dssp_output/")):
            subprocess.run("mkdssp -v " + ent_file + " " + dssp_file, shell=True)
        print(f"{pdb.upper()}'s dssp file has been generated.\n")
        with open(dssp_file, "r") as dssp:
            print(pdb)
            dssp = dssp.read()
            pattern = (r".....  ... . .  H....................................................................."
                       r"..................................................")
            for match in re.finditer(pattern, dssp):
                line = match.group()
                #print(line)
                helices.write(f"{line}\n")
                pattern = r".........."
                search = re.search(pattern,line)
                match = search.group()
                match = match.split(" ")[4]
                print(match)
                matches_list.append(match)
            while ("" in matches_list):
                matches_list.remove("")
            matches_list = np.array(matches_list).astype(int)
            for i in range(1,len(matches_list)):
                if matches_list[i] != int((matches_list[i-1])+1):
                    helices_list += 1
            #for index in range(0,len(helices_list)):
             #   helices_annotation.at[index, "pdb"] = pdb
              #  helices_annotation.at[index, "helix_number"] = index
                #helices_annotation.at[index, "helix_length"] = length
            helices.close()
            print(f"{pdb} has been annotated. Structure number: {counter}.")

    if option == "aleph":
        #try:
        structure = parser.get_structure(pdb, ent_file)
        for model in structure:
            for chain in model:
                for residue in chain:
                    atoms = residue.get_atoms()
                    atoms_list = {atom.get_name() for atom in atoms}
                    if all(name in atoms_list for name in ["N", "CA", "C", "O"]):
                        res_name = residue.get_resname()
                        residuesName_list.append(res_name)
                        segid = str(residue.get_segid)
                        regex = r"resseq=..."
                        match = re.search(regex, segid)
                        match = match.group()
                        match = match.split("=")[1]
                        match = match.split(" ")[0]
                        residuesPosition_list.append(match)
                        c = (str(chain).split('=')[1]).split('>')[0]
                        chain_list.append(c)
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
        helix_total, lengths = calculateHelicalContent(euclidean_distances, residuesName_list, residuesPosition_list, chain_list,pdb)
        print(pdb, len(lengths), helix_total)
        for helix in range(0,helix_total):
            helices_annotation.at[helix, "pdb"] = pdb
            length = lengths[helix]
            helices_annotation.at[helix, "helix_number"] = helix+1
            helices_annotation.at[helix, "helix_length"] = length
        print(f"{pdb} has been annotated. Structure number: {counter}.")
        helices_annotation.to_csv(f"{option}_{pdb}_helices_annotation.csv", sep="\t")
        helices_annotation.drop(helices_annotation.index, inplace=True)
        ca.clear(); x_ca.clear(); y_ca.clear(); z_ca.clear(); oxi.clear(); x_oxi.clear(); y_oxi.clear(); z_oxi.clear()
        euclidean_distances.clear(); residuesName_list.clear(); residuesPosition_list.clear(); chain_list.clear()
        #except:
        #print(f"Error: {pdb} PDB not found! Located in line {counter}.")

