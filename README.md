# Master-Thesis-Scripts

Hello!

This is the repository that I have created to store all the scripts developed during the elaboration of my Master Thesis Project.
All scripts have been separated in folders depending on their main usage.

Starting from the beginning, I downloaded the full PDB database using this script:

- <strong>download_fullPDB.sh</strong> (folder: fullPDB_filter): downloads each entry from the PDB database using the shell terminal. Later, gunzip of all .gz files is required to ensure all scripts correctly run.

Later, I started the filtering step of the data. This project has only used those structures identified as <strong>proteins</strong> and obtained through <strong>x-ray diffraction</strong> methods.
In order to obtain the desired data, I programmed the following script:

- <strong>filter_fullPDB.py</strong> (folder: full_PDB_filer): identifies and classifies entries depending on the experimental method and  molecular type labels.

After selecting only those entries of interest, I developed a script that allows for coiled-coil detection in structures. For that, I used SOCKET2 and DSSP softwares.

- <strong>socket_coiled_coils.py</strong> (folder: full_PDB_filter): creates DSSP files from PDB structures and, later, identifies coiled-coils with SOCKET2. All outputs from each software are stored in subfolders.

Additionally, I used the already existing code - developed previously by the team - that calculates characteristic vectors distances to predict secondary structure types, so I could annotate the helical, sheet and coil content from each structure.

- <strong>calculate_ss_content.py</strong>: uses characteristic vectors as mesure to identify secondary structure motifs. Further information can be found in the referenced literature at the end of this README).

Having filtered the PDB as I intented, I started working on the scripts for descriptors extraction. All the scripts intended to extract descriptors can be  located in the “descriptors” and “definitive” directories, and are further explained in my Master Thesis Project.

<strong>REFERENCES:</strong>

Team's previous work:

- Sammito, M., Millán, C., Rodríguez, D. D., De Ilarduya, I. M., Meindl, K., De Marino, I., Petrillo, G., Buey, R. M., De Pereda, J. M., Zeth, K., Sheldrick, G. M., & Usón, I. (2013). Exploiting tertiary structure through local folds for crystallographic phasing. Nature Methods 2013 10:11, 10(11), 1099–1101. https://doi.org/10.1038/nmeth.2644
- Medina, A., Triviño, J., Borges, R. J., Millán, C., Usón, I., & Sammito, M. D. (2020). ALEPH: a network-oriented approach for the generation of fragment-based libraries and for structure interpretation. Acta Crystallographica. Section D, Structural Biology, 76(Pt 3), 193. https://doi.org/10.1107/S2059798320001679

SOCKET2:

- Kumar, P., & Woolfson, D. N. (2021). Socket2: a program for locating, visualizing and analyzing coiled-coil interfaces in protein structures. Bioinformatics, 37(23), 4575. https://doi.org/10.1093/BIOINFORMATICS/BTAB631

Other:

- Empereur-Mot, C., Garcia-Seisdedos, H., Elad, N., Dey, S., & Levy, E. D. (2019). Geometric description of self-interaction potential in symmetric protein complexes. Scientific Data 2019 6:1, 6(1), 1–9. https://doi.org/10.1038/s41597-019-0058-x
