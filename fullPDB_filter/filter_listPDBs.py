#!/user/bin/env python
# Conda env: cctbx

# Usage: python script.py list(.txt) path_to_PDBs_directory
# This script filters PDBs into 3 outputs:
#   - Proteins obtained with x-ray diffraction.
#   - DNA/RNA molecules obtained with x-ray diffraction.
#   - Other PDBs that have been obtained with different methods.

# LIBRARIES IMPORTATION:
import sys
from Bio.PDB import *
parser = PDBParser(QUIET = True)
import iotbx.pdb.hierarchy
import pandas as pd

# ARGUMENTS:
pdb_list = list()
counter = 0
# Note: this path corresponds to where I have downloaded the full PDB. It must be changed to adapt to the user.
#fullPDB = os.getcwd()
list = sys.argv[1]
list = pd.read_csv(list, sep = '\t', header = 0)
fullPDB = sys.argv[2]
print(fullPDB)
results_file = open("XRAY_PROTEIN_PDBs.txt", "w")
notXRAY_file = open("notXRAY_PDBs.txt","w")
XRAY_DNA_RNA = open("XRAY_DNAorRNA_PDBs.txt","w")
obsolete = open("obsoletePDBs.txt","w")
unclassified = open("unclassified.txt","w")

# Iterate over files in input directory
for i in list.index:
    ele = list["pdb"][i]
    print(ele)
    try:
        ent_file = f"{fullPDB}pdb{ele}.ent"
        pdb = ele
        structure = parser.get_structure(pdb, ent_file)
        header = parse_pdb_header(ent_file)
        counter += 1
        is_protein = False
        is_na = False
        is_other = False
        chains_revised = []
        if "x-ray" in header["structure_method"]:
            pdb_in = iotbx.pdb.hierarchy.input(file_name=ent_file)
            for chain in pdb_in.hierarchy.only_model().chains():
                if chain.is_protein():
                    is_protein = True
                    chains_revised.append("protein")
                elif chain.is_na():
                    is_na = True
                    chains_revised.append("na")
            if "na" not in chains_revised:
                results_file.write(f"{pdb}\n")
            elif "na" in chains_revised:
                XRAY_DNA_RNA.write(f"{pdb}\n")
            else:
                unclassified.write(f"{pdb}\n")
        else:
            notXRAY_file.write(f"{pdb}\n")
        print(f"{pdb} has been classified. Structure number: {counter}.")
    except:
        counter+=1
        obsolete.write(f"{pdb}\n")
        print(f"Error: {pdb} PDB not found or obsolete! Located in line {counter}.")

results_file.close()
notXRAY_file.close()
XRAY_DNA_RNA.close()
obsolete.close()
unclassified.close()