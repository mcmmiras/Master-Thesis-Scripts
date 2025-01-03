# !usr/bin/env python

# Usage: python script.py fullPDB_path PDB_list(.txt)


# LIBRARIES IMPORTATION:
from Bio.SeqUtils import IsoelectricPoint
from Bio.PDB import *
parser = PDBParser(QUIET=True)
io = PDBIO()
sr = ShrakeRupley()
import os, sys, subprocess
import pandas as pd
import re
from Bio.SeqUtils import seq1
import traceback
import numpy as np
import math


# ARGUMENTS:
fullPDB = sys.argv[1]
list = sys.argv[2]
df = pd.read_csv(list, sep="\t")
df_results = pd.read_csv(list, sep="\t")
dir = os.getcwd()
if not "modified_PDBs" in os.listdir(dir):
    subprocess.run(f"mkdir modified_PDBs/", shell=True)
resultsHeptads = open("output_heptads_annotation.txt", "w")
errors = open(f"errors_heptads.txt", "w")
errors_list = []
array_generation = open("output_arrays_extraction.txt", "w")
data = pd.DataFrame(columns=["pdb","data","oligomer"])


# DICTIONARIES:
helicalStructure = {
    "hydrophobic": ["A", "V", "L", "I", "F", "W", "M", "C", "G"],  # Hydrophobic amino acids
    "hydrophilic": ["N", "Q", "S", "T"],  # Hydrophilic amino acids
    "charged": ["K", "R", "H", "D", "E"],  # Charged amino acids (acidic and basic)
    "polar": ["S", "T", "Q", "N"],  # Polar amino acids (often hydrophilic)
    "aromatic": ["F", "W", "Y"],  # Aromatic amino acids
    "disruptive": ["P"]  # Disruptive proline, it is hydrophobic but will be considered separatedly
}

aminoacids = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V", "MSE": "U"
}

oligomerization_states = {
	"monomeric": 1,
	"dimeric": 2,
    "trimeric":3,
    "tetrameric":4,
    "pentameric":5,
    "hexameric":6,
    "heptameric":7,
    "octameric":8,
    "nonameric":9
}

# FUNCTIONS:
# 1. Function to assign heptad repeats:
def assign_heptads(pdb, chains_seqs):
    heptads = str()
    for line in chains_seqs:
        i = 0
        #line = line.strip("X")
        if ">" not in line:
            seq = line
            print(f"Original sequence: {seq}")
            print(f"Original sequence length: {len(seq)}")
            seq = seq.strip("X")
            while i < (len(seq)-3):
                if seq[i] == "X":
                    i+=1
                    heptads = heptads+"-"
                else:
                    if (seq[i] in helicalStructure["hydrophobic"] and seq[i+3] in helicalStructure["hydrophobic"]):
                        heptads = heptads+"abcdefg"
                        i+=7
                    else:
                        i+=1
                        heptads = heptads+"-"
            heptads = heptads[:len(seq)]
            positions = [i.start() for i in re.finditer("abcdefg", heptads)]
            resultsHeptads.write(f"\n{pdb}\n{seq}\n{heptads}\n")
        # for ele in heptads:
        #     if ele == "a" or ele == "d":
        #         pos = heptads.index(ele)
        #         hydrophobicResidues = hydrophobicResidues+seq[pos]
            return seq, positions#, hydrophobicResidues
    #chains_seqs.close()


# 2. Function to change b-factors following heptad annotation:
# Saving results: Save the modified structure to a new PDB file:
def bfactorToHeptadAnnotation(structure, pdb, pos, seq):
    residuesList = []
    residuesElements = []
    abvResiduesList = str()
    foundHeptads = []
    for model in structure:
        for chain in model:
            if chain.id == "A":
                residues = chain.get_residues()
                for residue in residues:
                    residuesElements.append(residue)
                    if residue.get_resname() in "HOH":
                        residuesList.append(residue.get_resname())
    for position in pos:
        foundHeptads.append(seq[position:position+7])
    # for residue in residuesList:
    #     if residue in aminoacids:
    #         abvResiduesList = abvResiduesList + aminoacids[residue]
    #     else:
    #         abvResiduesList = abvResiduesList + "-"
    # # B-FACTORS MODIFICATION
    # residues_set = set()  # Track residues with modified B-factors
    # for heptad in foundHeptads:
    #     start_index = 0
    #     while (location := abvResiduesList.find(heptad, start_index)) != -1:
    #         for i in range(location, location + 7):
    #             residues_set.add(residuesElements[i])
    #             for atom in residuesElements[i]:
    #                 atom.set_bfactor(100 if i % 7 == 0 or i % 7 == 4 else 50)
    #         start_index = location + 1
    # for residue in residuesElements:
    #     if residue not in residues_set:
    #         for atom in residue:
    #             atom.set_bfactor(0)
    # io.set_structure(structure)
    # #if not f"modified_{pdb}.pdb" in os.listdir(f"{dir}/modified_PDBs/"):
    # io.save(f"{dir}/modified_PDBs/modified_{pdb}.pdb")
    return foundHeptads


# 3. Function to assign 'canonical heptad' score to each heptad:
# Reference: https://www.pnas.org/doi/10.1073/pnas.0604871103
# Canonical heptad:  Positions a and d are typically occupied by hydrophobic amino acids such as Leu, Ile, Val, and Ala,
# whereas residues at positions e and g are frequently polar or charged (5–8).
def assignCanonicalHeptadScore(heptads_list, positions):
    scores_list = []
    max_score = 0
    for ele in heptads_list:
        score = 0
        if ele[0] in helicalStructure["hydrophobic"]:
            score += 2
        if ele[4] in helicalStructure["hydrophobic"]:
            score += 2
        if ele[4] in helicalStructure["polar"] or ele[4] in helicalStructure["charged"]:
            score += 1
        if ele[6] in helicalStructure["polar"] or ele[6] in helicalStructure["charged"]:
            score += 1
        scores_list.append(score)
    for i, score in enumerate(scores_list):
        if score > max_score:
            max_score = score
            max_score_index = i
    best_heptad = heptads_list[max_score_index]
    best_heptad_pos = positions[max_score_index]
    return scores_list, best_heptad, max_score, best_heptad_pos


# 4. Function to assign 'environment hydrophobicity' of each residue in a 9-residue area:
def assignEnvironmentHydrophobicity(seq, initial_pos):
    scores_list = []
    i = initial_pos
    while i < len(seq):
        analyzed_positions = [i-4, i-3, i-1, i, i+1, i+3, i+4]
        score = 0
        max_analized_position = i+4
        for position in analyzed_positions:
            if position >= 0 and max_analized_position < len(seq):
                if seq[position] in helicalStructure["hydrophobic"]:
                    score += 1
                elif seq[position] in helicalStructure["polar"] or seq[position] in helicalStructure["charged"]:
                    score -= 1
                else:
                    score += 0
        scores_list.append(score)
        i += 1
    return scores_list


# EXECUTION:
# Iteration through all the PDBs list:
counter = 0
errors_num = 0
arrays_lists = []
for i in df.index:
    counter += 1
    try: # Obtain PDB from list
        pdb = df['pdb'][i]
        ent_file = f"pdb{pdb}.ent"
        ent_file = os.path.join(fullPDB, ent_file)
        try: # Assign heptad repeats
            structure = parser.get_structure(pdb,ent_file)
            for model in structure:
                for chain in model:
                    if Chain.Chain.get_id(chain) == "A":
                        chains_seqs = {chain.id:seq1(''.join(residue.resname for residue in chain))}
            sequence, positions = assign_heptads(pdb, chains_seqs.values())
            heptads = bfactorToHeptadAnnotation(structure, pdb, positions, sequence)
            canonicalScores, best_heptad, max_score, best_heptad_pos = assignCanonicalHeptadScore(heptads, positions)
            scores_hydrophobicity = assignEnvironmentHydrophobicity(sequence, best_heptad_pos)
            df.at[i, 'heptads'] = str(heptads)
            df.at[i, "heptads_pos"] = str(positions)
            df.at[i, "heptads_score"] = str(canonicalScores)
            df.at[i, "best_heptad"] = str(best_heptad)
            df.at[i, "best_heptad_pos"] = str(best_heptad_pos)
            df.at[i, "best_heptad_score"] = str(max_score)
            df_results.at[i, "initial_heptad"] = str(best_heptad)
            df_results.at[i, "initial_pos"] = str(best_heptad_pos)
            sequence_array = []
            scores_hydrophobicity_array = []
            print(f"Whole sequence without unknown residues: {sequence}")
            total_heptads = math.floor((len(sequence[best_heptad_pos:])/7))
            for i in range(0, total_heptads*7):
                sequence_array.append(sequence[best_heptad_pos:][i])
                scores_hydrophobicity_array.append(scores_hydrophobicity[i])
            sequence_array = np.array(sequence_array)
            scores_hydrophobicity_array = np.array(scores_hydrophobicity_array)
            df_results.at[i, "sequence_array"] = str(sequence_array)
            df_results.at[i, "hydrophobicity_array"] = str(scores_hydrophobicity_array)
            arrays_lists.append(scores_hydrophobicity_array)
            array_generation.write(
                f"\n{pdb.upper()}\n"
                f"Whole sequence without unknown residues: {sequence}\n"
                f"Most canonical-like heptad is '{best_heptad}' starting at position {best_heptad_pos+1} of the sequence.\n"
                f"Total found heptads: {total_heptads}.\n"
                f"Array of analyzed residues:\n{sequence_array}\n"
                f"Array of calculated hydrophobicity values:\n{scores_hydrophobicity_array}\n"
                f"Analyzed sequence array length: {len(sequence_array)}\n"
                f"Calculated values array length: {len(scores_hydrophobicity_array)}\n"
            )
            print(f"Most canonical-like heptad is '{best_heptad}' starting at position {best_heptad_pos+1} of the sequence.")
            print(f"Total found heptads: {total_heptads}.")
            print(f"Array of analyzed residues:\n{sequence_array}")
            print(f"Array of calculated hydrophobicity values:\n{scores_hydrophobicity_array}")
            print(f"Analyzed sequence array length: {len(sequence_array)}\n"
                  f"Calculated values array length: {len(scores_hydrophobicity_array)}")
            print(f"{pdb.upper()} has been correctly processed. Structure number {counter}.\n")
        except:
            print(f"Error in {pdb}: impossible to assign heptad repeats. Structure number {counter}.")
            errors_num += 1
            errors.write(f"{pdb}\n")
            errors_list.append(pdb)
            traceback.print_exc()
    except:
        errors_num += 1
        errors.write(f"{pdb}\n")
        errors_list.append(pdb)
        traceback.print_exc()
        print(f"Error in {pdb.upper()}. Could not process structure number {counter}.")

# CONVERT ALL DATA TO SAME LENGTH INPUTS:
arrays_lists = [array.tolist() for array in arrays_lists]
arrays_lists_length = [len(array) for array in arrays_lists]
arrays_mean = round(sum(arrays_lists_length) / len(arrays_lists_length))
arrays_mean = 28
print(arrays_mean)
#while arrays_mean%7 != 0:
    #arrays_mean += 1
    #print(arrays_mean)
print(arrays_mean)
larger = 0
for i in arrays_lists_length:
    if i > arrays_mean:
        diff = int((i - arrays_mean)/7)
        print(f"added rows: {diff}")
        larger += diff
normalization = np.zeros([len(arrays_lists_length)+larger, arrays_mean])
print(normalization)
print(normalization.shape)
row_displacement = 0
"""
for row, array in enumerate(arrays_lists):
    print(row)
    if len(array) == arrays_mean:
        normalization[int(row + row_displacement), :] = array
    elif len(array) < arrays_mean:
        padding = (arrays_mean - len(array)) // 2
        for col, ele in enumerate(array):
            normalization[int(row + row_displacement), int(col + padding)] = ele
    elif len(array) > arrays_mean:
        diff = math.ceil(len(array) / arrays_mean)
        for ind in range(diff):
            iter = arrays_mean * ind
            for col, ele in enumerate(array[iter:iter + arrays_mean]):
                normalization[int(row + row_displacement + ind), col] = ele
        row_displacement += diff - 1
"""
index_pdbs = 0
for row, array in enumerate(arrays_lists):
    print(len(array))
    pdb = df['pdb'][index_pdbs]
    if pdb in errors_list:
        index_pdbs += 1
    ent_file = os.path.join(fullPDB, f"pdb{df['pdb'][index_pdbs]}.ent")
    try:
        result = subprocess.check_output(f"grep 'AUTHOR DETERMINED BIOLOGICAL UNIT' {ent_file}", shell=True)
        result = str(result)
        for state in oligomerization_states.keys():
            if state.upper() in result:
                oligomer = oligomerization_states[state]
    except:
        try:
            result = subprocess.check_output(f"grep 'SOFTWARE DETERMINED QUATERNARY STRUCTURE:' {ent_file}", shell=True)
            result = str(result)
            for state in oligomerization_states.keys():
                if state.upper() in result:
                    oligomer = oligomerization_states[state]
        except:
            errors_num += 1
            errors.write(f"{pdb}\n")
            errors_list.append(pdb)
            traceback.print_exc()
            print(f"Error in {pdb.upper()}. Could not process structure number {counter}.")

    if len(array) == arrays_mean:
        for col, ele in enumerate(array):
            normalization[int(row+row_displacement)][col] = ele
        data.at[int(row+row_displacement), "pdb"] = df['pdb'][index_pdbs]
        data.at[int(row+row_displacement), "oligomer"] = oligomer
    if len(array) < arrays_mean:
        diff = (arrays_mean - len(array))/2
        if diff != int(diff):
            displacement = math.floor(diff)
            displacementF = displacement+1
        else:
            displacement = diff
            displacementF = diff
        for col, ele in enumerate(array):
            normalization[int(row+row_displacement)][int(col+displacement)] = ele
        data.at[int(row+row_displacement), "pdb"] = df['pdb'][index_pdbs]
        data.at[int(row + row_displacement), "oligomer"] = oligomer
    if len(array) > arrays_mean:
        diff = int((len(array) - arrays_mean)/7)
        print(f"Diff: {diff}")
        for ind in range(diff):
            i = 7*ind
            print(f"i value for iteration: {i}")
            normalization[int(row + row_displacement + ind), :arrays_mean] = array[i:i + arrays_mean]
            data.at[int(row + row_displacement + ind), "pdb"] = df['pdb'][index_pdbs]
            data.at[int(row + row_displacement + ind), "oligomer"] = oligomer
        remainder_start = 7 * diff
        normalization[int(row + row_displacement + diff), :arrays_mean] = array[
                                                                          remainder_start:remainder_start + arrays_mean]
        data.at[int(row + row_displacement + diff), "pdb"] = df['pdb'][index_pdbs]
        data.at[int(row + row_displacement + diff), "oligomer"] = oligomer
        row_displacement += diff
    index_pdbs += 1
data['data'] = normalization.tolist()
data.to_csv('data.csv', sep = '\t')
"""
        diff = math.ceil(len(array) / arrays_mean)
        for ind in range(0, diff):
            iter = arrays_mean*ind
            if len(array[array[iter:iter + arrays_mean]]) < arrays_mean:
                diff = len((array[array[iter:iter + arrays_mean]]) - len(array)) / 2
                if diff != int(diff):
                    displacement = math.floor(diff)
                    displacementF = displacement + 1
                else:
                    displacement = diff
                    displacementF = diff
                for col, ele in enumerate(array[array[iter:iter + arrays_mean]]):
                    normalization[int(row + row_displacement)][int(col + displacement)] = ele
            else:
                for col, ele in enumerate(array[iter:iter + arrays_mean]):
                    if col < arrays_mean:
                        normalization[int(row + row_displacement+ind)][col] = array[int(col+iter)]
"""

print(normalization)
print(normalization.shape)
print(larger)
print(f"\nOUTPUT SUMMARY:")
print(f"Structures correctly processed: {counter - errors_num}.")
print(f"Structures not correctly classified: {errors_num}.\n")

df.to_csv("output_descriptors_heptads.csv", sep="\t")
df_results.to_csv("output_descriptors.csv", sep="\t")
resultsHeptads.close()
errors.close()
array_generation.close()
