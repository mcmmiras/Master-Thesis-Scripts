# !/usr/bin/env python
import os
# Usage: python script.py pdb(name) pdb(file) hydrophobicity.scale(s)
# Available scales (the scale name must be introduced as described):	(---- recommended for alpha-helices)
	# ----kyte_doolittle: for identifying both hydrophobic surface-exposed regions as well as transmembrane regions. Regions with a positive value are hydrophobic.
	# hopp_woods: for identification of potentially antigenic sites in proteins. Apolar residues have been assigned negative values.
	# ----cornette: optimal hydrophobicity scale based on 28 published scales. Suitable for prediction of alpha-helices in proteins.
	# eisenberg: normalized consensus hydrophobicity scale which shares many features with the other hydrophobicity scales
	# rose: correlated to the average area of buried amino acids in globular proteins. Scale which is not showing the helices of a protein, but rather the surface accessibility.
	# janin: information about the accessible and buried amino acid residues of globular proteins.
	# ----engelman_ges: useful for predicting transmembrane regions in proteins. Similar to kyte_doolittle.
	# all: combination of all the available scales


# LIBRARIES:
# Importation of the required libraries:
import sys
from Bio.PDB import *
parser = PDBParser()
import matplotlib.pyplot as plt
import random
import os
import subprocess
import pandas as pd


# DICTIONARIES (HYDROPHOBICITY SCALES):
# We create dictionaries to store the hydropathy indexes of the residues based on different hydrophobicity scales:
kyte_doolittle = {
	"ALA": 1.8,
	"CYS": 2.5,
	"ASP": -3.5,
	"GLU": -3.5,
	"PHE": 2.8,
	"GLY": -0.4,
	"HIS": -3.2,
	"ILE": 4.5,
	"LYS": -3.9,
	"LEU": 3.8,
	"MET": 1.9,
	"ASN": -3.5,
	"PRO": -1.6,
	"GLN": -3.5,
	"ARG": -4.5,
	"SER": -0.8,
	"THR": -0.7,
	"VAL": 4.2,
	"TRP": -0.9,
	"TYR": -1.3
}

hopp_woods = {
	"ALA": -0.5,
	"CYS": -1.0,
	"ASP": 3.0,
	"GLU": 3.0,
	"PHE": -2.5,
	"GLY": 0.0,
	"HIS": -0.5,
	"ILE": -1.8,
	"LYS": 3.0,
	"LEU": -1.8,
	"MET": -1.3,
	"ASN": 0.2,
	"PRO": 0.0,
	"GLN": 0.2,
	"ARG": 3.0,
	"SER": 0.3,
	"THR": -0.4,
	"VAL": -1.5,
	"TRP": -3.4,
	"TYR": -2.3,
}

cornette = {
	"ALA": 0.2,
	"CYS": 4.1,
	"ASP": -3.1,
	"GLU": -1.8,
	"PHE": 4.4,
	"GLY": 0.0,
	"HIS": 0.5,
	"ILE": 4.8,
	"LYS": -3.1,
	"LEU": 5.7,
	"MET": 4.2,
	"ASN": -0.5,
	"PRO": -2.2,
	"GLN": -2.8,
	"ARG": 1.4,
	"SER": -0.5,
	"THR": -1.9,
	"VAL": 4.7,
	"TRP": 1.0,
	"TYR": 3.2
}

eisenberg = {
	"ALA": 0.62,
	"CYS": 0.29,
	"ASP": -0.90,
	"GLU": -0.74,
	"PHE": 1.19,
	"GLY": 0.48,
	"HIS": -0.40,
	"ILE": 1.38,
	"LYS": -1.50,
	"LEU": 1.06,
	"MET": 0.64,
	"ASN": -0.78,
	"PRO": 0.12,
	"GLN": -0.85,
	"ARG": -2.53,
	"SER": -0.18,
	"THR": -0.05,
	"VAL": 1.08,
	"TRP": 0.81,
	"TYR": 0.26
}

rose = {
	"ALA": 0.74,
	"CYS": 0.91,
	"ASP": 0.62,
	"GLU": 0.62,
	"PHE": 0.88,
	"GLY": 0.72,
	"HIS": 0.78,
	"ILE": 0.88,
	"LYS": 0.52,
	"LEU": 0.85,
	"MET": 0.85,
	"ASN": 0.63,
	"PRO": 0.64,
	"GLN": 0.62,
	"ARG": 0.64,
	"SER": 0.66,
	"THR": 0.70,
	"VAL": 0.86,
	"TRP": 0.85,
	"TYR": 0.76
}

janin = {
	"ALA": 0.3,
	"CYS": 0.9,
	"ASP": -0.6,
	"GLU": -0.7,
	"PHE": 0.5,
	"GLY": 0.3,
	"HIS": -0.1,
	"ILE": 0.7,
	"LYS": -1.8,
	"LEU": 0.5,
	"MET": 0.4,
	"ASN": -0.5,
	"PRO": -0.3,
	"GLN": -0.7,
	"ARG": -1.4,
	"SER": -0.1,
	"THR": -0.2,
	"VAL": 0.6,
	"TRP": 0.3,
	"TYR": -0.4
}

engelman_ges = {
	"ALA": 1.6,
	"CYS": 2.0,
	"ASP": -9.2,
	"GLU": -8.2,
	"PHE": 3.7,
	"GLY": 1.0,
	"HIS": -3.0,
	"ILE": 3.1,
	"LYS": -8.8,
	"LEU": 2.8,
	"MET": 3.4,
	"ASN": -4.8,
	"PRO": -0.2,
	"GLN": -4.1,
	"ARG": -12.3,
	"SER": 0.6,
	"THR": 1.2,
	"VAL": 2.6,
	"TRP": 1.9,
	"TYR": -0.7
}
total_scales = ["kyte_doolittle","hopp_woods","cornette","eisenberg","rose","janin","engelman_ges"]


# FUNCTIONS:
# Function for converting hydrophobicity scales into a 0-100 scale:
def hydrophobicity_converter(scale):
	# Variables with the maximum and minimum hydropathy indexes:
	if scale == rose:
		min_index, max_index = 0.52, 0.91
	else:
		min_index, max_index = 0, 0
		for i in scale:
			if scale[i] > max_index:
				max_index = scale[i]
			else:
				if scale[i] < min_index:
					min_index = scale[i]
	# Iteration through the dictionary to calculate and add the normalized values to the corresponding keys:
	for i in scale:
		# Normalization of the hydropathy indexes into a 0-100 scale:
		norm_value = ((scale[i]-min_index)/(max_index - min_index))*100
		if scale == hopp_woods:
			norm_value = (norm_value-100)*(-1)
		norm_value = round(norm_value, 1)
		# Transform values into a list, and append the new normalized value:
		scale[i] = [scale[i]]
		scale[i].append(norm_value)
	return scale

# Function for bfactors conversion into hydrophobic indexes:
def bfactors_converter(pdb_file, converted_scale):
# Iteration through the structure:
	for model in pdb_file:
		for chain in model:
			for residue in chain:
				# Obtention of the residue name abbreviated into 3-letter code:
				residue_3letters = Residue.Residue.get_resname(residue)
				# Iteration through the dictionary to find the key that matches the residue name:
				for x in converted_scale:
					if x == residue_3letters:
						# Assignment of the normalized hydropathy index to the bfactor variable:
						bfactor = converted_scale[x][1]
						# Substitution of all residue's atoms bfactor value with the new bfactor:
						for atom in residue:
							Atom.Atom.set_bfactor(atom, bfactor)


# Function to obtain the coordinates for later graphical representation:
def graphic_generation(pdb_file, ax, scale_name):
	# GRAPHICAL REPRESENTATION. Dictionary generation with the data to plot:
	# Generation of lists to save the residues' modified b-factors to obtain hydropathic indexes:
	sum_ind = 0
	x,y = [],[] # x is number of residue, y is bfactor
	i = 1
	for model in pdb_file:
		for chain in model:
			for residue in chain:
				for atom in residue:
					atom_name = Atom.Atom.get_name(atom)
					if atom_name=="CA":
						bfactor = Atom.Atom.get_bfactor(atom)
						# To store the data corresponding to each structure into their respective keys in the library
						x.append(i)
						y.append(bfactor) # Storing the b-factor values
						sum_ind += bfactor
						i+=1
	# HYDROPATHY SCORE calculation:
	hydropathy_score = round(sum_ind / len(x), 2)
	# GRAPHICAL REPRESENTATION. Plot generation with the stored data:
	# Plotting the stored data:
	ax.plot(x, y, color="#" + ''.join([random.choice('0123456789ABCDEF') for j in range(6)]), label=f"{sys.argv[2]}-{scale_name}")
	ax.axhline(50, linestyle="--", color="purple")
	# GRAPHICAL REPRESENTATION. Labelling and text in the plot:
	ax.legend(loc="lower right", title="Structure", ncol=1, fontsize="small", title_fontsize="medium")
	ax.set_title(f"Hydropathy profile: {scale_name} scale\n(hydropathy score: {hydropathy_score})", fontsize="11")
	ax.set_xlabel("Residues")
	ax.set_ylabel(f"Hydropathic index")
	ax.set_ylim(0,100)
	return hydropathy_score

def all_converter(pdb_file):
	fig, axes = plt.subplots(nrows=4, ncols=2, figsize=(12, 10))
	for i, scale in enumerate(total_scales):
		scale_name=total_scales[i]
		# Conversion to a 0-100 scale:
		converted_scale = hydrophobicity_converter(eval(scale_name))
		# bfactors conversion into hydrophobic indexes:
		bfactors_converter(pdb_file, converted_scale)
		# SAVING RESULTS. Save the modified structure to a new PDB file:
		io = PDBIO()
		io.set_structure(pdb_file)
		io.save(f"modified_{pdb_name}_{scale_name}.pdb")
		# Obtention of the coordinates for later graphical representation:
		ax = axes[i // 2, i % 2]
		graphic_generation(pdb_file, ax, scale_name)
	# To show the plot in the screen and store it:
	plt.delaxes(axes[3, 1])
	plt.tight_layout()  # Adjust subplot layout
	plt.savefig(f"hydropathy_{pdb_name}_multiple_scales.png")
	plt.show()


# ARGUMENTS:
# Hydrophobicity scale that will be applied to the pdb(file):
hydrophobicity_scale = sys.argv[3]
# PDB file taken and parsed:
pdb_files = sys.argv[2]
pdb_files = pd.read_csv(pdb_files, sep="\t")
fullPDB = sys.argv[1]
if not "output" in os.listdir(os.getcwd()):
	subprocess.run("mkdir output/", shell=True)
file = open("hydrophobicity_scores.csv","w")
file.write(f"pdb\thydrophobicity_score\n")

# CONVERSIONS:
if hydrophobicity_scale == 'all':
	# Conversion to a 0-100 scale and bfactors conversion into hydrophobic indexes for all available scales_
	for i in pdb_files.index:
		pdb_name = pdb_files.loc[i, "pdb"]
		print(pdb_name)
		ent_file = os.path.join(fullPDB, f"pdb{pdb_name}.ent")
		print(ent_file)
		pdb_file = parser.get_structure(pdb_name, ent_file)
		all_converter(pdb_file)
else:
	# Conversion to a 0-100 scale:
	converted_scale = hydrophobicity_converter(eval(hydrophobicity_scale))
	for i in pdb_files.index:
		pdb_name = pdb_files.loc[i, "pdb"]
		print(pdb_name)
		ent_file = os.path.join(fullPDB, f"pdb{pdb_name}.ent")
		print(ent_file)
		pdb_file = parser.get_structure(pdb_name, ent_file)
		# bfactors conversion into hydrophobic indexes:
		bfactors_converter(pdb_file, converted_scale)
		# SAVING RESULTS:
		# Save the modified structure to a new PDB file:
		io = PDBIO()
		io.set_structure(pdb_file)
		io.save(f"output/modified_{pdb_name}_{hydrophobicity_scale}.pdb")
		# Obtention of the coordinates for graphical representation:
		fig, ax = plt.subplots()
		hydrophobicity_score = graphic_generation(pdb_file, ax, hydrophobicity_scale)
		plt.savefig(f"output/hydropathy_{pdb_name}_{hydrophobicity_scale}.png")
		#plt.show()
		file.write(f"{pdb_name}\t{hydrophobicity_score}\n")


print("Your plot with the hydropathy profile of your structure has been generated successfully.")