#!/user/bin/env python
# Conda env: cctbx

# Usage: python script.py path_to_PDBs_directory
# This script filters PDBs into 3 outputs:
#   - Proteins obtained with x-ray diffraction.
#   - DNA/RNA molecules obtained with x-ray diffraction.
#   - Other PDBs that have been obtained with different methods.

# LIBRARIES IMPORTATION:
import os
import sys
from Bio.PDB import *
parser = PDBParser(QUIET = True)
import iotbx.pdb.hierarchy

# ARGUMENTS:
pdb_list = list()
counter = 0
fullPDB = sys.argv[1]
# OUTPUTS:
results_file = open("XRAY_PROTEIN_PDBs.txt", "w")
notXRAY_file = open("notXRAY_PDBs.txt","w")
XRAY_DNA_RNA = open("XRAY_DNAorRNA_PDBs.txt","w")
obsolete = open("obsoletePDBs.txt","w")
unclassified = open("unclassified.txt","w")

# Iterate over files in input directory: everything is stored into txt files.
for ele in os.listdir(fullPDB):
    if ".ent" in ele:
        try:
            # Parse PDB (.ent file):
            ent_file = os.path.join(fullPDB,ele)
            pdb = ele[3:7]
            structure = parser.get_structure(pdb, ent_file)
            header = parse_pdb_header(ent_file)
            counter += 1
            is_protein = False
            is_na = False
            is_other = False
            chains_revised = list()
            # Firstly, a separation between experimental method: we select those obtained with x-ray diffraction:
            if "x-ray" in header["structure_method"]:
                pdb_in = iotbx.pdb.hierarchy.input(file_name=ent_file)
                # Then, classification of each structure into protein, nucleic acid or other:
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