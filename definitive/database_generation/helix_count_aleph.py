#!/user/bin/env python

# Usage: python script.py path_to_fullPDB working_directory helical_percent(decimal)

import os
import sys
import re
import pandas as pd
import subprocess
import traceback
from Bio.PDB import *
parser = PDBParser(QUIET=True)

# LISTS:
aminoacids = [
    "ALA",
    "ARG",
    "ASN",
    "ASP",
    "CYS",
    "GLN",
    "GLU",
    "GLY",
    "HIS",
    "ILE",
    "LEU",
    "LYS",
    "MET",
    "PHE",
    "PRO",
    "SER",
    "THR",
    "TRP",
    "TYR",
    "VAL",
    "MSE"
]


# ARGUMENTS:
fullPDB = sys.argv[1]
annotations_directory = sys.argv[2]
helical_percent = sys.argv[3]
dir = os.getcwd()
if "output_summary" not in os.listdir(annotations_directory):
    subprocess.run(f"mkdir output_summary", shell=True)
errors = open("output_summary/errors_pdb.csv", "w")
processed = open("output_summary/processed_pdb.csv","w")
results = pd.DataFrame(columns=["pdb","helix_count","helix_residues","total_residues","helix_residues/total_residues"])
discarded = pd.DataFrame(columns=["pdb","helix_count","helix_residues","total_residues","helix_residues/total_residues"])

# EXECUTION:
counter = 0
errors_num = 0
for index, pdb_directory in enumerate(os.listdir(annotations_directory)):
    if len(pdb_directory) == 4:
        counter += 1
        helix_total = 0
        residues_helix_total = 0
        total_residues = 0
        try:
            ent_file = os.path.join(fullPDB,f"pdb{pdb_directory}.ent")
            structure = parser.get_structure(f"{pdb_directory}", f"{ent_file}")
            residues = structure.get_residues()
            for residue in residues:
                if residue.get_resname() in aminoacids:
                    total_residues += 1
            with open(f"{annotations_directory}/{pdb_directory}/pdb{pdb_directory}_input_search.pdb") as pdb_file:
                helix_extraction = open(f"{annotations_directory}/{pdb_directory}/pdb{pdb_directory}_helix_annotation.txt","w")
                pdb_lines = pdb_file.read()
                pattern = r"^.*\bHELIX\b.*$"
                matches = re.finditer(pattern, pdb_lines, re.MULTILINE)
                for match in matches:
                    helix_extraction.write(f"{match.group()}\n")
                    helix_total += 1
                    pattern2 = r"                      ..\d$"
                    match2 = re.search(pattern2, match.group())
                    match2 = match2.group().replace(" ", "")
                    residues_helix_total += int(match2)
            pdb_file.close()
            helix_extraction.close()
            processed.write(f"{pdb_directory}\n")
            helical_content = round(int(residues_helix_total)/int(total_residues),2)
            if helical_content >= float(helical_percent):
                results.at[index,"pdb"] = pdb_directory
                results.at[index,"helix_count"] = int(helix_total)
                results.at[index,"helix_residues"] = int(residues_helix_total)
                results.at[index, "total_residues"] = int(total_residues)
                results.at[index, "helix_residues/total_residues"] = helical_content
            else:
                discarded.at[index, "pdb"] = pdb_directory
                discarded.at[index, "helix_count"] = int(helix_total)
                discarded.at[index, "helix_residues"] = int(residues_helix_total)
                discarded.at[index, "total_residues"] = int(total_residues)
                discarded.at[index, "helix_residues/total_residues"] = helical_content
            print(f"Annotation of {pdb_directory.upper()} file has been made correctly. Structure number: {counter}.")
        except:
            errors_num += 1
            errors.write(f"{pdb_directory}\n")
            traceback.print_exc()
            print(f"Error in {pdb_directory.upper()}: annotation was not possible. Structure number: {counter}.")

processed.close()
errors.close()
results.to_csv(f"output_summary/descriptor_helixcount{helical_percent}_aleph.csv",sep="\t")
discarded.to_csv(f"output_summary/discarded_helixcount{helical_percent}_aleph.csv",sep="\t")

print(f"\n")
print(f"OUTPUT SUMMARY:")
print(f"Total correctly annotated structures: {counter-errors_num}.")
print(f"Total failed annotations: {errors_num}.")