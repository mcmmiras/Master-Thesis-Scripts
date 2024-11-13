#!/usr/bin/env python
import os
from unittest.mock import inplace

# Usage: pymol -qc script.py
# Usage: python script.py fullPDB_path PDB_list(.txt)

#   -r: opens pymol GUI
#   -qc: prints the results in the terminal but does not open the pymol GUI

# LIBRARIES:
from Bio.PDB import *
parser = PDBParser(QUIET=True)
import pandas as pd
from math import sqrt, log
from pymol import *
import sys
import subprocess
import matplotlib.pyplot as plt

# ARGUMENTS:
fullPDB = sys.argv[1]
list = sys.argv[2]
df_list = pd.read_csv(list, sep="\t")
dir = os.getcwd()
if not "output_rASA_files" in os.listdir(dir):
    subprocess.run(f"mkdir output_rASA_files/", shell=True)
if not "output_split_chains" in os.listdir(dir):
    subprocess.run(f"mkdir output_split_chains/", shell=True)
if not "output_stickiness_scales" in os.listdir(dir):
    subprocess.run(f"mkdir output_stickiness_scales/", shell=True)
errors = []


# FUNCTIONS:
# Function 1: create dataframe for each PDB.
def create_dataframe():
    df_pdb = pd.DataFrame(data={"residue_name": [], 'residue_number': [], 'chain': [], 'rASA_complexed': [], 'rASA_unbound': [],
                                "delta": [], ">25_rASAc": [], ">25_rASAu": [], "location": []})
    df_pdb["residue_name"] = df_pdb["residue_name"].astype(str)
    df_pdb["residue_number"] = df_pdb["residue_number"].astype(str)
    df_pdb["chain"] = df_pdb["chain"].astype(str)
    df_pdb[">25_rASAc"] = df_pdb[">25_rASAc"].astype(str)
    df_pdb[">25_rASAu"] = df_pdb[">25_rASAu"].astype(str)
    df_pdb["location"] = df_pdb["location"].astype(str)
    return df_pdb

# Function 2: perform PyMol chains extraction and rASA calculation per residue.
def extract_chains_rASA(ent_file, df_pdb,pdb,residues_list):
    cmd.load(ent_file)  # Loading the PDB file
    obj = f"pdb{pdb}"
    rASA_complexed = cmd.get_sasa_relative(obj, quiet=0, _self=cmd)  # Calculating rASA in complexed state
    # Iterate through all rASA values and store them:
    i = 0
    for key in rASA_complexed.keys():
        chain = key[2]
        residue_num = key[3]
        rASA = round(rASA_complexed[key], 2)
        df_pdb.at[i, "residue_name"] = residues_list[i]
        df_pdb.at[i, "residue_number"] = residue_num
        df_pdb.at[i, "chain"] = chain
        df_pdb.at[i, "rASA_complexed"] = rASA
        if rASA > 0.25:
            df_pdb.at[i, ">25_rASAc"] = "yes"
        else:
            df_pdb.at[i, ">25_rASAc"] = "no"
        i += 1
    cmd.split_chains()  # Split all chains into individual objects
    chains_id = cmd.get_chains()  # Get the chains identifiers
    # Iterate through all split chains and store the new rASA values:
    i = 0
    for chain in chains_id:
        obj = f"pdb{pdb}_{chain}"
        rASA_unbound = cmd.get_sasa_relative(obj, quiet=0, _self=cmd)  # Calculating rASA in unbound state per chain
        cmd.save(f"./output_split_chains/pdb{pdb}_{chain}.pdb", f"pdb{pdb}_{chain}")  # Output with the split chain
        # Iteration through all calculated rASA values for unbound states:
        for key in rASA_unbound.keys():
            rASA = round(rASA_unbound[key], 2)
            df_pdb.at[i, "rASA_unbound"] = rASA
            if rASA > 0.25:
                df_pdb.at[i, ">25_rASAu"] = "yes"
            else:
                df_pdb.at[i, ">25_rASAu"] = "no"
            i += 1
    cmd.delete("all")
    return df_pdb

# Function 3: calculate delta value of residues' rASA between complexed and unbound states of a protein.
def calculate_delta():
    for i, row in df_pdb.iterrows():
        delta = sqrt((df_pdb.loc[i, "rASA_unbound"] - df_pdb.loc[i, "rASA_complexed"])**2)
        df_pdb.at[i, "delta"] = round(delta,1)

# Function 4: calculate delta value of residues' rASA between complexed and unbound states of a protein.
def assign_location():
    for i, row in df_pdb.iterrows():
        if df_pdb.loc[i, ">25_rASAc"] == "yes" and df_pdb.loc[i, "delta"] == 0:
            df_pdb.at[i, "location"] = "surface"
        if df_pdb.loc[i, ">25_rASAc"] == "no" and df_pdb.loc[i, "delta"] > 0 and df_pdb.loc[i, ">25_rASAu"] == "yes":
            df_pdb.at[i, "location"] = "interface"

#Function 5: calculate log of location propensity per residue.
def calculate_residues_stickiness(df_pdb, pdb, oligomer):
    # First, we will create the dataframe to store the stickiness scale:
    df_res = pd.DataFrame(data={"residue": ["ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE", "LYS",
                                                         "LEU", "MET", "ASN", "PRO", "GLN", "ARG", "SER", "THR", "VAL",
                                                         "TRP",
                                                         "TYR"]})
    df_res.set_index('residue', inplace=True)
    df_res["freq_interface"] = 0
    df_res["freq_surface"] = 0
    df_res["log(freq_interface/freq_surface)"] = 0
    df_res["log_normalized"] = 0
    # Then, we will identify the PDB file and read its corresponding output with rASA values:
    df_pdb = pd.read_csv(os.path.join(os.getcwd(),f"output_rASA_files/{df_pdb}"), sep="\t", header=(0))
    # Iteration to calculate frequencies and log between interface/surface:
    for i, row in df_pdb.iterrows():
        residue = df_pdb.iloc[i, 1]
        location = df_pdb.iloc[i,9]
        if location == "interface":
            df_res.at[residue,"freq_interface"] = df_res.loc[residue,"freq_interface"]+1
        elif location == "surface":
            df_res.at[residue,"freq_surface"] = df_res.loc[residue,"freq_surface"]+1
    # Values to later generate a normalizing scale:
    min_log = 0
    max_log = 0
    for i,row in df_res.iterrows():
        if df_res.loc[i,"freq_interface"] != 0 and df_res.loc[i,"freq_surface"] !=0:
            log = (math.log(df_res.loc[i,"freq_interface"]/df_res.loc[i,"freq_surface"]))
            df_res.at[i,"log(freq_interface/freq_surface)"] = log
        else:
            log = 0
        # Values to later generate a normalizing scale:
        if log > max_log:
            max_log = log
        elif log < min_log:
            min_log = log
    # Transformation of values into a -1,0,+1 scale range:
    for i, row in df_res.iterrows():
        log = df_res.loc[i,"log(freq_interface/freq_surface)"]
        log_norm = (log - min_log) / (max_log - min_log)
        log_norm = 2 * log_norm - 1
        df_res.at[i,"log_normalized"] = log_norm
    # Graphical representation of the calculated stickiness scale:
    plt.figure(figsize=(10, 10))
    df_res.plot(y="log_normalized", use_index=True, kind="bar")
    plt.axhline(0, color='k')
    plt.title(f"{pdb}: {oligomer}", loc="center")
    plt.savefig(f"./output_stickiness_scales/{pdb}_stickiness_scale.png")
    # Finally, we will save our results:
    df_res.to_csv(f"./output_stickiness_scales/{pdb}_stickiness_scale.csv", sep = "\t")



# EXECUTION:
pymol.finish_launching(["pymol","-qc"]) # Run PyMol from script without opening the GUI
# Iterate all PDBs from a given list:
for i in df_list.index:
    pdb = df_list['pdb'][i]
    try:
        if f"{pdb}_rASA.csv" not in os.listdir(os.path.join(os.getcwd(),"output_rASA_files")):
            df_pdb = create_dataframe()
            ent_file = f"pdb{pdb}.ent" # General use for fullPDB
            ent_file = os.path.join(fullPDB, ent_file)
            residues_names = []
            # Parsing of PDB file for storing the residue names:
            try:
                structure = parser.get_structure(pdb, ent_file)
                residues_obj = Structure.Structure.get_residues(structure)
                for res in residues_obj:
                    if res != "HOH":
                        residues_names.append(Residue.Residue.get_resname(res))
            except:
                print(f"Error in parsing {pdb}.")
            try:
                # Extraction of all chains separatedly and calculation of rASA in complexed and unbound states:
                df_pdb = extract_chains_rASA(ent_file,df_pdb,pdb,residues_names)
            except:
                print(f"Error in {pdb}: could not extract chains.")
                errors.append(f"{pdb}\n")
            try:
                # Calculate delta (differences) between rASAc and rASAu values:
                calculate_delta()
            except:
                print(f"Error in {pdb}: could not calculate delta.")
                errors.append(f"{pdb}\n")
            try:
                # Assign residues' location within the folded structure:
                assign_location()
            except:
                print(f"Error in {pdb}: could not assign location.")
                errors.append(f"{pdb}\n")
            # Store results in separated csv files for each PDB:
            df_pdb.to_csv(f"./output_rASA_files/{pdb}_rASA.csv", sep="\t")
    except:
        print(f"Error in {pdb}: could not create dataframe.")
        errors.append(f"{pdb}\n")

# Terminate PyMol:
cmd.quit()

# Return PDB that have generated errors:
with open("./unprocessed_PDBs.txt","w") as file:
    file.writelines(errors)
    for error in errors:
        error = error.split(f"\n")[0]
        subprocess.run(f"rm -r ./output_rASA_files/{error}_rASA.csv", shell=True)
    file.close()

# Iterate through rASA output for each PDB file:
df_list.set_index("pdb", inplace=True)
for df_pdb in os.listdir(os.path.join(os.getcwd(), f"output_rASA_files")):
    try:
        pdb = df_pdb[:4]
        oligomer = df_list.loc[pdb,"oligomer"]
        if not f"{pdb}_stickiness_scale.csv" in os.listdir(os.path.join(os.getcwd(), f"output_stickiness_scales")):
            # Calculate the stickiness scale for each PDB file:
            calculate_residues_stickiness(df_pdb, pdb, oligomer)
    except:
        print(f"Error: the {df_pdb} file was not found.")