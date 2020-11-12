#!/usr/bin/python

try:
    fasta_seqs = open("/home/nmoreyra/Documents/Repositorios/GeneticDiversityEstimation/data/dna_example.fasta")
except IOError:
	print("File does not exit!")

seqs = {}
name = ''
for line in fasta_seqs:
#	print(line)
	line = line.rstrip() #let's discard the newline at the end (if any)
	#now we need to distinguish header from sequence
	if line[0] == '>': #identify the id of each sequence
		words = line.split()
        	name = words[0][1:]
        	seqs[name] = ''
	else:
	        seqs[name] = seqs[name] + line

nucl = {'N':'N' 'A':'A', 'T':'T', 'C':'C', 'G':'G'}

ws = 10 #window size
s = 1 #step to slide the window
for name,seq in seqs.items():
	slen = len(seq) #sequence length
	pos = 0 #iteration variable to walk over the sequence
	nw = int(slen/ws)+1 #number of windows in the total sequence length
	while slen-pos-s >= ws-1: #while the windows does not exceed the end of the sequence
		pC = seq[pos:pos+ws].count('C')/slen
		pG = seq[pos:pos+ws].count('G')/slen
		pT = seq[pos:pos+ws].count('T')/slen
		pA = seq[pos:pos+ws].count('A')/slen
print()no_c=dna_sequence.count('c')
pos=pos+s

