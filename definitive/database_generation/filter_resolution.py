#!/user/bin/env python

# Usage: python script.py path_to_PDBs_directory listofPDBs
# This scripts filters structures with resolution no higher than 4 Ångström.

# LIBRARIES IMPORTATION:
from Bio.PDB import parse_pdb_header
import os, sys
import pandas as pd
import traceback

# ARGUMENTS:
fullPDB = sys.argv[1]
pdb_list = sys.argv[2]
pdb_list = pd.read_csv(pdb_list, sep="\t", header=None)
pdb_list = pdb_list.set_index(pdb_list[0])
selected = open("selected_pdb.csv", "w")
selected.write(f"pdb\tresolution\n")
discarded = open("discarded_pdb.csv", "w")
discarded.write(f"pdb\tresolution\n")
no_resolution = open("no_resolution_pdb.csv", "w")
no_resolution.write(f"pdb\tobtention_method\n")
errors = open("errors_pdb.csv", "w")


# EXECUTION:
for pdb in pdb_list.index:
    ent_file = os.path.join(fullPDB,f"pdb{pdb}.ent")
    try:
        result = parse_pdb_header(ent_file)
        resolution = result['resolution']
        method = result["structure_method"]
        try:
            if not resolution:
                no_resolution.write(f"{pdb}\t{method}\n")
            else:
                if resolution <= 4.0:
                    selected.write(f"{pdb}\t{resolution}\n")
                else:
                    discarded.write(f"{pdb}\t{resolution}\n")
            print(f"{pdb.upper()} was correctly classified.")
        except:
            print(f"Error in {pdb}: classification was not possible.")
            errors.write(f"{pdb}\terror type: classification\n")
            traceback.print_exc(file=errors)
            errors.write(f"\n")

    except:
        print(f"Error in {pdb}: could not parse PDB header.")
        errors.write(f"{pdb}\terror type: parsing\n")
        traceback.print_exc(file=errors)
        errors.write(f"\n")

selected.close()
discarded.close()
no_resolution.close()
errors.close()