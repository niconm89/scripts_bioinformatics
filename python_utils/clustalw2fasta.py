#!/bin/python

#%% Imports
import argparse
from Bio import SeqIO

def convert_aln2fasta(infile, outfile):
	records = SeqIO.parse(infile, "clustal")
	dict_records = {}
	for rec in records:
		dict_records[rec.id] = rec.seq

	with open(outfile, 'wt') as OUT:
		for record in sorted(dict_records.keys()):
			OUT.write(">"+record+"\n"+str(dict_records[record])+"\n")

def main():
	parser = argparse.ArgumentParser(
		description='''clustalw2fasta.py conversts a clustalw alignment file (.aln) into a fasta file.''',
		epilog="""End of the help""")

	parser.add_argument('-f', '--infile', type=str, required=True, help='Full path to the clustalw alignment file to convert')
	parser.add_argument('-o', '--outfile', type=str, required=True, help='Full path where the fasta file will be placed')

	args=parser.parse_args()

	convert_aln2fasta(args.infile,args.outfile)

#%% Main program
if __name__ == '__main__':
    main()
