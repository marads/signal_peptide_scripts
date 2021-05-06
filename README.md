# signal_peptide_scripts

# check_Fo_SP.py
# python 3
# Author: Mara de Sain 
# Date: June 2019
# input: fasta file with signal peptide as header, e.g. the file all_putative_effectors_concatenated_clustered.fasta from Peters FoEC.py script can be used

# Python script to check if a signal peptide identified by SignalP adheres to the following three rules created by Martijn Rep, if so it is considered a 'True' signal peptide and the gene will be used for further analysis:
# 1. Size: between 15 and 30 amino acids
# 2. Processing site: the last three residues are small-any-small
# 3. Hydrophobicity: contains a stretch of 9 to 20 consecutive hydrophobic residues.
# This stretch:
# a. starts at least after the 1st and at most after the 8th residue
# b. ends at most 2 residues before the last three residues
# c. contains no hydrophilic residues
# d. contains a maximum of 3 hydroxylated residues

# Summary:
# [M][1-7any][9-20hydrophobicwithmax3hydroxylated][0-2any]-small-any-small

# any = GALMFWKQESPVICYHRNDT
# hydrophobic = AVILMFWCGP
# hydrophilic = RKDEQNH
# hydroxylated = STY
# small = GAVSTC
