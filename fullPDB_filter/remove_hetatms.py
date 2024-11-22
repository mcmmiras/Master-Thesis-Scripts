#!/user/bin/env python

# Usage: python script.py fullPDB_path PDB_list option[filter/remove]

import os
import sys
from urllib.parse import non_hierarchical

import pandas as pd
import subprocess
import traceback

# ARGUMENTS:
pdb_list = sys.argv[2]
pdb_list = pd.read_csv(pdb_list, sep="\t", header=None)
pdb_list = pdb_list.set_index(pdb_list[0])
print(pdb_list)
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
for pdb in pdb_list.index:
    ent_file = os.path.join(fullPDB,f"pdb{pdb}.ent")
    if option == "filter":
        try:
            found = subprocess.check_output(f"grep -q HETATM {ent_file}", shell=True)
            discarded.write(f"{pdb}\n")
            print(f"{pdb.upper()} has been discarded.")
        except:
            selected.write(f"{pdb}\n")
            print(f"{pdb.upper()} has been selected.")
            error.write(f"{pdb}\n")
            traceback.print_exc(file=error)
            error.write(f"\n")
    if option == "remove":
        try:
            found = subprocess.check_output(f"grep -q HETATM {ent_file}", shell=True)
            subprocess.run(f"grep -v HETATM {ent_file} > 'modified_PDBs/pdb{pdb}.ent'", shell=True)
            modified.write(f"{pdb}\n")
            print(f"HETATMs from {pdb.upper()} have been removed correctly.")
        except:
            non_modified.write(f"{pdb}\n")
            print(f"HETATMs from {pdb.upper()} were not found.")
            error.write(f"{pdb}\n")
            traceback.print_exc(file=error)
            error.write(f"\n")

error.close()
non_modified.close()
modified.close()
discarded.close()
selected.close()