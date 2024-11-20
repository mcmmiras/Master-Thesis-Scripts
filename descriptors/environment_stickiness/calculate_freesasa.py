#!/usr/bin/env python

# Usage: pymol -qc script.py
# Usage: python script.py fullPDB_path PDB_list(.txt)
# Conda env: freesasa

#   -r: opens pymol GUI
#   -qc: prints the results in the terminal but does not open the pymol GUI

# LIBRARIES:
#from Bio.PDB import *
#parser = PDBParser(QUIET=True)
import pandas as pd
from math import sqrt, log
from Bio.PDB import *

parser = PDBParser(QUIET=True)
from pymol import *
import freesasa
freesasa.Parameters.defaultParameters["algorithm"] = "ShrakeRupley"
print(f"Freesasa software default parameters: {freesasa.Parameters.defaultParameters}")
import sys, os
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

# Function 2: perform PyMol chains extraction.
def extract_chains(ent_file, pdb):
    split_chains = []
    cmd.load(ent_file)  # Loading the PDB file
    cmd.split_chains()  # Split all chains into individual objects.
    chains_id = cmd.get_chains()  # Get the chains identifiers
    if not f"{pdb}" in os.listdir(os.path.join(dir,"output_split_chains")):
        subprocess.run(f"mkdir ./output_split_chains/{pdb}/", shell=True)
    for chain in chains_id:
        cmd.save(f"./output_split_chains/{pdb}/pdb{pdb}_{chain}.pdb", f"pdb{pdb}_{chain}")  # Output with the split chain
        split_chains.append(f"./output_split_chains/{pdb}/pdb{pdb}_{chain}.pdb")
    cmd.delete("all")
    return split_chains

# Function 3a and 3b: perform rASA calculation per residue in complexed (3a) and unbound (3b) states.
def calculate_rASA_complexed(ent_file, df_pdb, residues_names):
    structure = freesasa.Structure(ent_file)
    result = freesasa.calc(structure)
    rASA = result.residueAreas()
    i = 0
    for chain in rASA.keys():
        for residue in rASA[chain]:
            #print(chain, residue, round(rASA[chain][residue].total, 2))
            value = rASA[chain][residue].total
            df_pdb.at[i, "residue_name"] = residues_names[i]
            df_pdb.at[i, "residue_number"] = residue
            df_pdb.at[i, "chain"] = chain
            df_pdb.at[i, "rASA_complexed"] = round(value, 2)
            if value > 0.25:
                df_pdb.at[i, ">25_rASAc"] = "yes"
            else:
                df_pdb.at[i, ">25_rASAc"] = "no"
            i+=1

def calculate_rASA_unbound(split_chains, df_pdb):
    i = 0
    for pdb in split_chains:
        structure = freesasa.Structure(pdb)
        result = freesasa.calc(structure)
        rASA = result.residueAreas()
        for chain in rASA.keys():
            for residue in rASA[chain]:
                # print(chain, residue, round(rASA[chain][residue].total, 2))
                value = rASA[chain][residue].total
                df_pdb.at[i, "residue_number"] = residue
                df_pdb.at[i, "chain"] = chain
                df_pdb.at[i, "rASA_unbound"] = round(value, 2)
                if value > 0.25:
                    df_pdb.at[i, ">25_rASAu"] = "yes"
                else:
                    df_pdb.at[i, ">25_rASAu"] = "no"
                i += 1

# Function 4: calculate delta value of residues' rASA between complexed and unbound states of a protein.
def calculate_delta(df_pdb):
    for i, row in df_pdb.iterrows():
        delta = sqrt((df_pdb.loc[i, "rASA_unbound"] - df_pdb.loc[i, "rASA_complexed"])**2)
        df_pdb.at[i, "delta"] = round(delta,1)

# Function 5: calculate delta value of residues' rASA between complexed and unbound states of a protein.
def assign_location(df_pdb):
    for i, row in df_pdb.iterrows():
        if df_pdb.loc[i, ">25_rASAc"] == "yes" and df_pdb.loc[i, "delta"] == 0:
            df_pdb.at[i, "location"] = "surface"
        if df_pdb.loc[i, ">25_rASAc"] == "no" and df_pdb.loc[i, "delta"] > 0 and df_pdb.loc[i, ">25_rASAu"] == "yes":
            df_pdb.at[i, "location"] = "interface"

# Function 6: calculate log of location propensity per residue.
# First, we will create the dataframe to store the stickiness scale:
def calculate_residues_stickiness(df_pdb, pdb, oligomer):
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
    min_log = -1
    max_log = +1
    for i,row in df_res.iterrows():
        if df_res.loc[i,"freq_interface"] != 0 and df_res.loc[i,"freq_surface"] !=0:
            log = (math.log(df_res.loc[i,"freq_interface"]/df_res.loc[i,"freq_surface"]))
            df_res.at[i,"log(freq_interface/freq_surface)"] = log
            print(df_res.loc[i,"log(freq_interface/freq_surface)"])
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
        df_res.at[i,"log_normalized"] = round(log_norm,3)
    # Graphical representation of the calculated stickiness scale:
    plt.figure(figsize=(10, 10))
    df_res.plot(y="log_normalized", use_index=True, kind="bar")
    plt.axhline(0, color='k')
    plt.title(f"{pdb}: {oligomer}", loc="center")
    plt.savefig(f"./output_stickiness_scales/{pdb}_stickiness_scale.png")

# EXECUTION:
pymol.finish_launching(["pymol","-qc"]) # Run PyMol from script without opening the GUI
# Iterate all PDBs from a given list:
for i in df_list.index:
    pdb = df_list['pdb'][i]
    print(pdb)
    residues_names = []
    try:
        df_pdb = create_dataframe()
        ent_file = f"pdb{pdb}.ent" # General use for fullPDB
        ent_file = os.path.join(fullPDB, ent_file)
        try:
            structure = parser.get_structure(pdb, ent_file)
            residues_obj = Structure.Structure.get_residues(structure)
            for res in residues_obj:
                if Residue.Residue.get_resname(res) != "HOH":
                    residues_names.append(Residue.Residue.get_resname(res))
        except:
            print(f"Error in parsing {pdb}.")
        try:
            # Calculation of each individual residue rASA:
            calculate_rASA_complexed(ent_file, df_pdb, residues_names)
            # Separate each chain individually:
            split_chains = extract_chains(ent_file, pdb)
            # Calculation of each individual residue rASA, now per each single chain:
            calculate_rASA_unbound(split_chains, df_pdb)
        except:
            print(f"Error in {pdb}: could not calculate rASA.")
    except:
        print(f"Error in {pdb}: could not create dataframe.")
        errors.append(f"{pdb}\n")
    try:
        # Calculate delta (differences) between rASAc and rASAu values:
        calculate_delta(df_pdb)
    except:
        print(f"Error in {pdb}: could not calculate delta.")
        errors.append(f"{pdb}\n")
    try:
        # Assign residues' location within the folded structure:
        assign_location(df_pdb)
    except:
        print(f"Error in {pdb}: could not assign location.")
        errors.append(f"{pdb}\n")
    # Store results in separated csv files for each PDB:
    df_pdb.to_csv(f"./output_rASA_files/rASA_{pdb}.csv", sep="\t")

# Return PDB that have generated errors:
with open("./unprocessed_PDBs.txt","w") as file:
    file.writelines(errors)
    for error in errors:
        error = error.split(f"\n")[0]
        subprocess.run(f"rm -r ./output_rASA_files/output_{error}.csv", shell=True)
    file.close()


# Terminate PyMol:
cmd.quit()

# Iterate through rASA output for each PDB file:
df_list.set_index("pdb", inplace=True)
for df_pdb in os.listdir(os.path.join(os.getcwd(), f"output_rASA_files")):
    print(df_pdb)
    try:
        pdb = df_pdb[5:9]
        #oligomer = df_list.loc[pdb,"oligomer"]
        #print(oligomer)
        oligomer = "" # For cases in which no previous oligomer classification has been made
        if not f"{pdb}_stickiness_scale.csv" in os.listdir(os.path.join(os.getcwd(), f"output_stickiness_scales")):
            # Calculate the stickiness scale for each PDB file:
            calculate_residues_stickiness(df_pdb, pdb, oligomer)
            print("Finished")
    except:
        print(f"Error: the {df_pdb} file was not found.")
