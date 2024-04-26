#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script: gbk2fasta.py
Author: Nicolás Nahuel Moreyra (niconm89@gmail.com)
GitHub: @niconm89
Date: Mon Aug  1 10:06:21 2022

Description:
    
Usage:
    python gbk2fasta.py

Requirements:
    - Biopython, time
"""

#%% Imports
import argparse
from Bio import SeqIO
from time import time

#%% Functions definition
def revcomp(sequence):
	complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
	reverse_complement = "".join(complement.get(base, base) for base in reversed(sequence))
	return reverse_complement
#end
def dbk2fasta(GBK, FASTA):
	with open(GBK, 'rt') as input_handle, open(FASTA, 'wt') as output_handle:
		sequences = SeqIO.parse(input_handle, "genbank")
		count = SeqIO.write(sequences, output_handle, "fasta")
		print("Converted %i records" % count)
#end
def dbk2fasta_split(GBK, FASTA):
	dbk2fasta(GBK, FASTA)
	mt_sequence = ''
	with open(FASTA, 'rt') as INPUT:
		for seq_record in SeqIO.parse(INPUT, "fasta"):
			mt_sequence = str(seq_record.seq)

	coordinates = {}
	with open(GBK, 'rt') as input_handle, open(FASTA, 'wt') as output_handle:
		for seq_record in SeqIO.parse(input_handle, "genbank"):
			for seq_feature in seq_record.features:
				if seq_feature.type == "gene":
					if 'product' in seq_feature.qualifiers:
						name = seq_feature.qualifiers['product'][0]
					else:
						name = seq_feature.qualifiers['gene'][0]
					location = str(seq_feature.location).split(":")
					start = int(location[0][1:])
					end = int(location[1].split("]")[0])
					strand = str(location[1].split("]")[1][1])
					coordinates[name] = [start, end, strand]

		#print(coordinates)
		for feature, coords in coordinates.items():
			#coordinates[name] = [start, end, strand]
			output_handle.write(">" + feature + "\n")
			if coords[2] == '-':
				reverse_gene = revcomp(mt_sequence[coords[0]:coords[1]])
				output_handle.write(reverse_gene + "\n")
			else:
				output_handle.write(mt_sequence[coords[0]:coords[1]] + "\n")
#end
#%% Menu -> is executed when the script is called independently
def usage():
	parser = argparse.ArgumentParser(
		description='''gbk2fasta.py extract sequences from a gbk file into a multifasta file''',
		epilog="""End of the help""")

	parser.add_argument('-i', '--gbkfile', type=str, required=True, help='Path to the genbank input file.')
	parser.add_argument('-o', '--outfasta', type=str, required=True, help='Path to save output fasta file with sequences.')
	parser.add_argument('-s', '--split', action='store_true', required=False, default=False, help='Generate a multifasta file with sequences generated from annotations.')

	return parser.parse_args()
#end
#%% Main program
if __name__ == '__main__':
	args = usage()
	start = time()
	if args.split:
		dbk2fasta_split(args.gbkfile, args.outfasta)
	else:
		dbk2fasta(args.gbkfile, args.outfasta)
	print(f'Time taken to run: {time() - start} seconds.')
#%% End
