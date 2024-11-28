#!/user/bin/env python

# Usage: python script.py fullPDB_path PDB_list option[filter/remove]

import os
import sys
import re
import pandas as pd
import subprocess
import traceback

# ARGUMENTS:
pdb_list = sys.argv[2]
pdb_list = pd.read_csv(pdb_list, sep="\t", header=0)
fullPDB = sys.argv[1]
option = sys.argv[3]

# OUTPUTS:
# Option 1 [filter]: filter PDBs containing HETATMs (heteroatoms o non-standard atoms) and keep those not containing HETATMs
if option == "filter":
    discarded = open("filter_discardedPDBs.txt","w")
    selected = open("filter_selectedPDBs.txt","w")
    error = open("filter_errorPDBs.txt","w")
# Option 2 [remove]: remove all HETATMs lines from a PDB
if option == "remove":
    modified = open("remove_modifiedPDBs.txt", "w")
    non_modified = open("remove_non_modifiedPDBs.txt","w")
    error = open("remove_errorPDBs.txt","w")
    if "modified_PDBs/" not in os.listdir(os.getcwd()):
        subprocess.run(f"mkdir modified_PDBs/" , shell=True)

# EXECUTION:
for i in pdb_list.index:
    pdb = pdb_list["pdb"][i]
    ent_file = os.path.join(fullPDB,f"pdb{pdb}.ent")
    # if option == "filter":
    #     try:
    #         with open(f"{ent_file}", "r") as file:
    #             file = file.read().splitlines()
    #             for line in file:
    #                 pattern = r"HETATM"
    #                 match = re.search(pattern, line)
    #                 if str(match) != "None":
    #                     pattern = r"HETATM...........MSE"
    #                     match = re.search(pattern, line)
    #                     if str(match) != "None":
    #                         print(f"{pdb.upper()} has been selected.")
    #                     else:
    #                         if not pdb in discarded.read():
    #                             discarded.write(f"{pdb}\n")
    #                         print(f"{pdb.upper()} has been discarded.")
    #                 else:
    #                     if not pdb in selected.read():
    #                         selected.write(f"{pdb}\n")
    #     except:
    #         error.write(f"{pdb}\n")
    #         traceback.print_exc(file=error)
    #         error.write(f"\n")

    if option == "remove":
        try:
            with open(f"{ent_file}","r") as file:
                file = file.read().splitlines()
                new_ent = open(f"modified_PDBs/pdb{pdb}.ent","w")
                extracted = open(f"modified_PDBs/heteroatoms_{pdb}.ent", "w")
                for line in file:
                    pattern = r"HETATM"
                    match = re.search(pattern, line)
                    modified.write(f"{pdb}\n")
                    if str(match) != "None":
                        pattern = r"HETATM...........MSE"
                        match = re.search(pattern, line)
                        if str(match) != "None":
                            new_ent.write(f"{line}\n")
                        else:
                            extracted.write(f"{line}\n")
                    else:
                        new_ent.write(f"{line}\n")
            print(f"HETATMs from {pdb.upper()} have been removed correctly.")
        except:
            non_modified.write(f"{pdb}\n")
            print(f"HETATMs from {pdb.upper()} were not found.")
            error.write(f"{pdb}\n")
            traceback.print_exc(file=error)
            error.write(f"\n")

try:
    error.close()
    non_modified.close()
    modified.close()
    discarded.close()
    selected.close()
except:
    pass

'''
numberCC = subprocess.run(f"grep 'COILED COILS PRESENT' {out_file} | cut -c25-30", shell=True,
                                              capture_output=True, text=True)
'''