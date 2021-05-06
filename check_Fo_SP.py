# check_Fo_SP.py
# python 3
# Author: Mara de Sain 
# Date: June 2019
# input: fasta file with signal peptide as header, e.g. the file all_putative_effectors_concatenated_clustered.fasta from Peters FoEC.py script can be used

# Python script to check if a signal peptide identified by SignalP adheres to the following three rules created by Martijn Rep, if so it is considered a 'True' signal peptide:
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


import os, sys, re
from Bio import SeqIO
from argparse import ArgumentParser


def max_3_hydroxylated_residues(possible_hydrophobic_stretch):
	# checks rule 3d
	search_hydroxylated_aa = re.findall("[STY]", possible_hydrophobic_stretch)
	if len(search_hydroxylated_aa) < 4:
		return True
	else:
		return False


def no_hydrophilic_residues(possible_hydrophobic_stretch):
	# checks rule 3c
	search_hydrophilic_aa = re.findall("[RKDEQNH]", possible_hydrophobic_stretch)
	if len(search_hydrophilic_aa) == 0:
		return True
	else:
		return False


def hydrophobic_stretch_present(possible_SP):
	possible_hydrophobic_stretches = []
	if len(possible_SP)<21:
		possible_hydrophobic_stretches.append(possible_SP[-12:-3])
		possible_hydrophobic_stretches.append(possible_SP[-13:-4])
		possible_hydrophobic_stretches.append(possible_SP[-14:-5])
		for sequence in possible_hydrophobic_stretches:
			#checks rules 3c and 3d
			if no_hydrophilic_residues(sequence) == True and max_3_hydroxylated_residues(sequence) == True:
				return True
	elif len(possible_SP) == 21:
		possible_hydrophobic_stretches.append(possible_SP[-13:-3])
		possible_hydrophobic_stretches.append(possible_SP[-13:-4])
		possible_hydrophobic_stretches.append(possible_SP[-14:-5])
		for sequence in possible_hydrophobic_stretches:
			#checks rules 3c and 3d
			if no_hydrophilic_residues(sequence) == True and max_3_hydroxylated_residues(sequence) == True:
				return True
	elif len(possible_SP) == 22:
		possible_hydrophobic_stretches.append(possible_SP[-14:-3])
		possible_hydrophobic_stretches.append(possible_SP[-14:-4])
		possible_hydrophobic_stretches.append(possible_SP[-14:-5])
		for sequence in possible_hydrophobic_stretches:
			#checks rules 3c and 3d
			if no_hydrophilic_residues(sequence) == True and max_3_hydroxylated_residues(sequence) == True:
				return True
	else:
		extra_aa_hydrophobic_stretch = len(possible_SP) - 22
		start_1 = -3 - 9 - extra_aa_hydrophobic_stretch
		start_2 = -4 - 9 - extra_aa_hydrophobic_stretch
		start_3 = -5 - 9 - extra_aa_hydrophobic_stretch
		possible_hydrophobic_stretches.append(possible_SP[start_1:-3])
		possible_hydrophobic_stretches.append(possible_SP[start_2:-4])
		possible_hydrophobic_stretches.append(possible_SP[start_3:-5])
		for sequence in possible_hydrophobic_stretches:
			#checks rules 3c and 3d
			if no_hydrophilic_residues(sequence) == True and max_3_hydroxylated_residues(sequence) == True:
				return True


def sp_length_ok(possible_SP):
	# checks rule 1
	if len(possible_SP)>14 and len(possible_SP)<31:
		return True
	else:
		return False


def adheres_to_SP_rules(SP_seq):
	# checks rules 2 and 3a+b
	find_re_possible_SP = re.findall("[M][GALMFWKQESPVICYHRNDT]{0,7}[AVILMFWCGPSTY]{9,20}[GALMFWKQESPVICYHRNDT]{0,2}[GAVSTC][GALMFWKQESPVICYHRNDT][GAVSTC]", SP_seq)
	if len(find_re_possible_SP) == 1:
		possible_SP = find_re_possible_SP[0]
		return sp_length_ok(possible_SP) and hydrophobic_stretch_present(possible_SP)
	elif len(find_re_possible_SP) > 1:
		print("Something went wrong: multiple signal peptides were identified")
		return False


def is_input_aa(SP_seq):
	# check that all input characters are amino acids
	search_non_aa = re.search("[^GALMFWKQESPVICYHRNDT]", SP_seq)
	if search_non_aa != None:
		print("The signal peptide " + SP_seq + " contains a non amino acid character: " + search_non_aa.group())
		return False
	else:
		return True


def write_correct_SP_to_file(input_path, output_path):
	with open(output_path, "a+") as output_file, open(input_path, "r") as input_file:
		for record in SeqIO.parse(input_file, "fasta"):
			SP_seq = record.id.upper().strip()
			# if input contains non amino acid characters, no other functions are called by the following line
			if is_input_aa(SP_seq) == True and adheres_to_SP_rules(SP_seq) == True:
				SeqIO.write(record, output_file, "fasta")
				print("The signal peptide: " + SP_seq + " does adhere to Martijns rules, Yeah!")
			else:
				print("The signal peptide: " + SP_seq + " does not adhere to Martijns rules")


def get_input_path(args):
	if args.input_path == None:
		print("\n[ERROR!]\nPlease specify a fasta file containing signal peptides as headers using: -i [infile]")
		print("Use -h to see additional options.\n")
		quit()
	else:
		return args.input_path


def get_paths(args):	
	input_path = get_input_path(args)
	output_path = args.output_path
	return (input_path, output_path)


def init_argument_parser():
	usage = "\n" + "*** python correct_SP6b.py -i [fasta file with signal peptide sequences as headers] <options> ***" 
	description = "*** Description: Check if signal peptides adhere to the following rule: [M][1-7any][9-20hydrophobicwithmax3hydroxylated][0-2any]-small-any-small. If so, print sequence in fasta format to output file. ***"
	parser = ArgumentParser(description = description, usage=usage)
	parser.add_argument("-i", "--in", dest = "input_path", help = "provide the absolute path to a fasta file with putative signal peptide sequences as headers")
	parser.add_argument("-o", "--out", dest = "output_path", help = "absolute output_path; default = ~/Desktop/output.fasta", default = os.path.expanduser("~/Desktop/output.fasta"))
	return parser


def main():
	parser = init_argument_parser()
	args = parser.parse_args()
	(input_path, output_path) = get_paths(args)
	write_correct_SP_to_file(input_path, output_path)


main()

