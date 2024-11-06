# !usr/bin/env python

# Usage: python script.py fullPDB_path PDB_list(.txt)


# LIBRARIES IMPORTATION:
from Bio.PDB import *
parser = PDBParser(QUIET=True)
io = PDBIO()
#from biotite import *
import os, sys, subprocess
import pandas as pd
import re


# ARGUMENTS:
fullPDB = sys.argv[1]
list = sys.argv[2]
df = pd.read_csv(list, sep="\t")
dir = os.getcwd()
if not "output_dssp_files" in os.listdir(dir):
    subprocess.run(f"mkdir output_dssp_files/", shell=True)
if not "output_rASA_files" in os.listdir(dir):
    subprocess.run(f"mkdir output_rASA_files/", shell=True)
# if not "modified_PDBs" in os.listdir(dir):
#     subprocess.run(f"mkdir modified_PDBs/", shell=True)


# FUNCTIONS:
# Function 1: calculate surface and interface of given protein
def calculate_rASA(ent_file, pdb, structure):
    rASA = open(f"./output_rASA_files/output_rASA_{pdb}.txt", "w")
    rASA.write(f"chain\tresidue\trASA_complexed\trASA_unbound\n")
    dssp_file = f"./output_dssp_files/{pdb}.dssp"
    if not f"{pdb}.dssp" in os.listdir(f"{dir}/output_dssp_files/"):
        subprocess.run("mkdssp -v " + ent_file + " " + dssp_file, shell=True)  # Create dssp file
    dssp = DSSP(structure[0], dssp_file)
    print("yes")
    for key in dssp.keys():
        # key contains the residue identifier
        residue = dssp[key]
        chain = key[0]
        residue_name = residue[1]
        residue_sasa = round(residue[3],2)  # SASA is typically at index 3
        rASA.write(f"{chain}\t{residue_name}\t{residue_sasa}\n")
    rASA.close()

    #return biotite.structure.sasa(structure[0])
# EXECUTION:
# Iteration through all the PDBs list:
for i in df.index:
    # 0. OBTAIN PDB FROM LIST
    try:
        pdb = df['pdb'][i]
        #ent_file = f"pdb{pdb}.ent" # General use for fullPDB
        ent_file = f"{pdb}.pdb" # For testing
        ent_file = os.path.join(fullPDB, ent_file)
        structure = parser.get_structure(pdb, ent_file)
        # 1. STICKINESS SCALE OF THE PROTEIN
        #   1.1. Calculate surface and interface:
        try:
            calculate_rASA(ent_file, pdb, structure)
            #df.at[i, "net_charge"] = calculate_net_charge(ent_file)
        except:
           print(f"Error in {pdb}: cannot calculate surface area.")
    except:
       print(f"Error: {pdb} PDB not found!")


#df.to_csv("output_descriptors.csv", sep="\t")


"""
Need to program a function to extract all chains from complexed PDB and then save a new PDB for each chain.
Then, calculate the difference between residues in complexed/unbound state. Extract delta.
Automate everything.

"""