#!/usr/bin/python
#NUcleotide Diversity Estimation : NUDE

def gde(window_size,step, infile):
	from Bio import AlignIO

	infile = "/home/nmoreyra/Documents/Repositorios/GeneticDiversityEstimation/data/Drosophila_boankobu2.aln" # quitar luego
	alignment = AlignIO.read(infile,"clustal")

	seqs = {}
	ids = {}
	cont = 0
	for line in alignment:
		seqs[cont] = line.seq
		ids[cont] = line.id
		cont+=1

	pi = []
	ss = len(seqs) #dictionary size
	com_seq = (ss*(ss-1))/2
	ws = 100 # = window size
	st = 10 # = step #to slide the window
	lenseq = 0
	n=0
	mut_in_win = []
	for name,seq in seqs.items():
		if n==0:
			lenseq = len(seq)
			if (lenseq%st)!= 0:
				numwin = int(lenseq/st)
			else:
				numwin = (lenseq/st)-1
		for i in range(name+1,ss):
			next_seq = seqs[i]
			#space to think
			v = 0-st
			ext = ws
			nwindow = 0
			while v+ext <= lenseq:
				window1 = seq[v+st:v+ws]
				window2 = next_seq[v+st:v+ws]
				v += st
				#
				window_mutations = 0	
				for j in range(len(window2)):
					if window1[j] != window2[j]:
						window_mutations+=1
				if n > 0:
					mut_in_win[nwindow]+=window_mutations
					nwindow+=1
				else:
					mut_in_win.append(window_mutations)
	#			print(window1,window2,window_mutations)
			
				if v+ext+ext>lenseq:
					ext = v+ext+ext-lenseq
		n+=1
	
	for k in range(len(mut_in_win)):
		pi.append(mut_in_win[k]/com_seq/ws)			
	#print(mut_in_win)
	#print(pi)
	return(pi)















