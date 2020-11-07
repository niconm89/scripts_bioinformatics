#!/usr/bin/python

#%% Imports
import argparse
import os
import re


#%% Function definitions
def convert_fasta2phylip(infile, outfile, mode, aln_format):
	dict_seqs = {}
	with open(infile,'rt') as inputfile:
		FASTA = inputfile.readlines()
		if FASTA[0][0] != '>':
			print("Input file does not have fasta format.\nLeaving...")
			exit(1)
		for line in FASTA:
			line = line.rstrip() #let's discard the newline at the end (if any)
			#now we need to distinguish header from sequence
			if line[0] == '>': #identify the id of each sequence
				name = check_nameseq(mode,line[1:]) #avoid ">"
				dict_seqs[name] = ''
			else:
				dict_seqs[name] = dict_seqs[name] + line
	if check_len(dict_seqs):
		if check_alphabet(dict_seqs):
			name_seqs = list(dict_seqs.keys())
			nseqs = len(name_seqs)
			len_aln = len(dict_seqs[name_seqs[0]])
			header = do_header(nseqs, len_aln)
			do_alignment(dict_seqs, mode, aln_format, header, outfile)
			
		else:#sequences have not allowed characters
			print("Leaving...\n")
			exit(1)		
	else:#sequences have different lengths
		print("Leaving...\n")
		exit(1)
#end convert_fasta2phylip

def convert_phylip2fasta(infile, outfile, mode, nogaps = False):
	with open(outfile, 'wt') as OUT:
		seqs_dict = seqs_from_phylip(infile)
		for record in seqs_dict.items():
			if nogaps:
				OUT.write(">" + record[0] + "\n" + record[1].replace('-','') + "\n")
			else:
				OUT.write(">" + record[0] + "\n" + record[1] + "\n")

def seqs_from_phylip(infile):
	dict_seqs = {}
	with open(infile,'rt') as input_file:
		HEADER = input_file.readline() #read the header
		PHYLIP = input_file.readlines() #read the rest of the file
		header = HEADER.strip().split(' ')
		nseqs = int(header[0])
		nchars = int(header[1])
		nlines = sum([ 1 for line in PHYLIP if line.strip() != '' ])
		if PHYLIP[0].strip() == '': #control por empty line between header and body
			PHYLIP.pop(0)
		if nlines == nseqs: #sequential format
			for line in PHYLIP:
				line = line.strip()
				if line != '':
					taxon_id = line[:-nchar].rstrip()
					taxon_seq = line[-nchars:]
					dict_seqs[taxon_id] = taxon_seq
		else: #interleaved format
			seqs_order = {}
			i = 0
			while i < nseqs: #get the first part (with id) of each taxon
				i += 1
				line = PHYLIP.pop(0).strip() #remove this line from file
				seqs_order[i] = line
				dict_seqs[i] = ''
			#Read the rest of the lines with only sequences
			count = 0
			for line in PHYLIP: #get the rest of the lines for each taxon
				line = line.rstrip()
				if line == '':
					continue
				else:
					taxon_seq = line.replace(' ','')
					count += 1
					dict_seqs[count] += taxon_seq
					if count == nseqs:
						count = 0
			acumulate_len = len(dict_seqs[1])
			#back to the first lines
			rest_len = nchars - acumulate_len #number of characters in each first line
			for i in range(1,nseqs+1):
				first_part_taxon = seqs_order[i].replace(' ','')
				taxon_id = first_part_taxon[:-rest_len]
				dict_seqs[i] = first_part_taxon[-rest_len:] + dict_seqs[i]
				dict_seqs[taxon_id] = dict_seqs.pop(i) #remove key = i and update key taxon_id			
	return dict_seqs
#End seqs_from_phylip

def do_alignment(seqs_dict, format_mode, format_aln, header, outfile):
	with open(outfile, 'wt') as OUT:
		OUT.write(header+"\n") #write the header
		
		if format_aln == "sequential":
			for record in seqs_dict.items():
				seq = spaces_in_sequence(record[1])
				id_seq = record[0]
				id_len = len(id_seq)
				line_aln = ''
				if format_mode == 'strict':
					nspaces = 10 - id_len #at this position id_len <= 10
					line_aln = id_seq + ' '*nspaces + seq
				else:
					line_aln = id_seq + ' ' + seq
				OUT.write(line_aln+"\n")
		else: #alignment format = interleaved
			id_seqs = list(seqs_dict.keys())
			lenseq = int(header.split(' ')[1])
			max_idlen = max([ len(idseq) for idseq in id_seqs ]) #max ID len
			for idseq in id_seqs:
				id_len = len(idseq)
				seq = seqs_dict[idseq][0:20]
				seq_region = spaces_in_sequence(seq)
				firs_line_seq = ''
				if format_mode == 'strict':
					nspaces = 10 - id_len #at this position id_len <= 10
					first_line_seq = idseq + ' '*nspaces + seq_region
				else: #format_mode = 'relaxed'
					nspaces = max_idlen - id_len
					first_line_seq = idseq + ' ' + seq_region
				OUT.write(first_line_seq+"\n")		
			OUT.write("\n")			
			i = 0
			while i < lenseq-1:
				if i == 0:
					i = 20
				else:
					i += 30
				for idseq in id_seqs:
					seq = seqs_dict[idseq][i:i+30]
					seq_region = spaces_in_sequence(seq)
					OUT.write(seq_region+"\n")
				OUT.write("\n")			
#end do_alignment

def spaces_in_sequence(sequence):
	aux_seq = ''
	for i in range(0,len(sequence),10):
		aux_seq += sequence[i:i+10] + ' '
	return	aux_seq[:-1]	

def do_header(nseqs, nchars):
	return str(nseqs) + ' ' + str(nchars)

def check_alphabet(seqs_dict):
	checking = True
	alphabet = ['A','G','C','T','U','Y','R','W','S','K','M','B','D','H','V','X','N','?','O','-']
	for seq in seqs_dict.items():
		for i,character in enumerate(seq[1].upper()):
			if character not in alphabet:
				checking = False
				print("Character not recognized at position " + str(i+1) + "\n")
				break
	return checking	

def check_len(seqs_dict):
	correct = True
	true_length = 0
	for i,record in enumerate(seqs_dict.items()):
		if i == 0:
			true_length = len(record[1])
		else:
			if len(record[1]) != true_length:
				correct = False
				print("Sequences must have the same length.")
				break
	return correct		

def check_nameseq(mode,id_seq):
	if '\n' in id_seq:
		print("Newline character not allowed in sequence ID.\nLeaving...")
		exit(1)
	if mode == "strict":
		if len(id_seq) > 10:
			id_seq = id_seq[0:10]
	#if mode = relaxed, id_seq will have no modification
	return id_seq

def check_files(infile,outfile):
	#infile
	if not os.path.exists(infile):
		print("Error: Input file does not exists! " + infile)
		exit(1)
	#outfile
	if outfile != 'convertion.phylip':
		#cwd = os.getcwd() + '/' #Get the current working directory and add slash (/)
		if '/' in outfile:
			elements = outfile.split('/')
			path_to_dir = '/'.join(elements[:-1])
			if not os.path.isdir(path_to_dir):
				print("Error: Output folder does not exist! " + path_to_dir)
				exit(1)

def usage():
	parser = argparse.ArgumentParser(
		description='''fasta2phylip.py converts a fasta alignment file into a phylip format.\nUser can select strict or relaxed format.''',
		epilog="""For more information, read the Phylip format documentation at https://evolution.genetics.washington.edu/phylip/doc/sequence.html""")

	parser.add_argument('--infile', metavar='<Input File>', type=str, required=True, help='Path to the fasta alignment file to convert. All sequences must have the same length (number of characters).')
	parser.add_argument('--outfile', metavar='<Output File>', type=str, required=False, default = 'convertion.phylip', help='Path where the phylip format file will be placed. By default, the output file will be placed in the current working directory. [default: convertion.phylip]')
	parser.add_argument('--mode', type=str, required=False, choices=['strict', 'relaxed'], default='strict', help='Convertion mode. In strict mode, id sequence have to be shorter than 10, otherwise it will be truncated. [default: strict]')
	parser.add_argument('--format', type=str, required=False, choices=['sequential', 'interleaved'], default='sequential', help='Alignment format. [default: sequential]')
	parser.add_argument('--reverse', action='store_true', required=False, help='phylip2fasta mode. The program will expect a phylip format file to convert to a fasta file. If not indicated (--mode), the input file is assumed to have a strict phylip format. If taxon name has a space, e.g. H. sapiens, it will be remove like H.sapiens.')
	parser.add_argument('--nogaps', action='store_true', required=False, help='Filter out gaps from the alignment. If declare, is mandatory to state the reverse mode (--reverse).')

	args=parser.parse_args()
	
	return args

def main():
	args = usage()
	check_files(args.infile,args.outfile)
	if not args.reverse:
		if args.nogaps:
			print("--nogaps option only compatible with reverse mode.")
			exit(1)		
		convert_fasta2phylip(args.infile,args.outfile,args.mode,args.format)
	else:
		if args.nogaps:
			convert_phylip2fasta(args.infile,args.outfile,args.mode, True)
		else:
			convert_phylip2fasta(args.infile,args.outfile,args.mode)
#end of main function

#%% Main program
if __name__ == '__main__':
	main()
	

