#!/home/nmoreyra/Soft/miniconda3/bin/python

#%% Imports
import os
import sys
#import gzip
from Bio import SeqIO

#%%Functions
def dict_pep(pep_path,biopython=False):
        'pep_path is the location to the predicted pepides file (generated with transdecoder).'
        #print("Creating peptides dictionary.\n")
        seqs = {}
        with open(pep_path,"rt") as pepfile:
                if biopython == False:
                        peptidome = pepfile.readlines()
                        for line in peptidome:
                                line = line.rstrip() #let's discard the newline at the end (if any)
                                if len(line) > 0 and line[0] == '>': #identify the id of each sequence
                                        words = line.split()
                                        name = words[0][1:]
                                        seqs[name] = ''
                                else:
                                        seqs[name] = seqs[name] + line
                else:
                        PEP = SeqIO.parse(pepfile,'fasta')
                        for i,pep in enumerate(PEP):
                                seqs[str(pep.id)] = str(pep.seq)
        return seqs
# end
def get_genes_multiples_prots(GFF_FILE):
	'Receive a GFF file from NCBI and return an ID list of genes with the number proteins annotated because of different isoforms.'
	count_genes = {}
	gene2prot = {}
	with open(GFF_FILE,'rt') as GFF:
		LINES = GFF.readlines()
		for line in LINES:
			columns = line.split("\t")
			if columns[2] == "gene":
				col9_list = columns[8].split(";")
				protid = col9_list[0][3:]
				for part in col9_list:
					if "gene=" in part:
						geneid = part.split("=")[1]
						if geneid not in count_genes:
							count_genes[geneid] = 1
							gene2prot[geneid] = [protid]
						else:
							count_genes[geneid] += 1
							gene2prot[geneid].append(protid)
	return count_genes, gene2prot
#end
def find_longest_prot(PROT_LIST, PROTSEQS):
	maxlen = 0
	maxid = ''
	for prot in PROT_LIST:
		len_prot = len(PROTSEQS[prot])
		if len_prot > maxlen:
			maxlen = len_prot
			maxid = prot
	return maxid
#end
######################################################################################################
#%% Main program

arguments = sys.argv

gff = arguments[1]
print(gff.split(".")[0])
#gff = "/home/nmoreyra/Data/cbuzz_genomes/09_detox-genes/drosophila_datasets/gff_faa_files/Ccap.gff"
prots = arguments[2]
#prots = "/home/nmoreyra/Data/cbuzz_genomes/09_detox-genes/drosophila_datasets/gff_faa_files/Ccap.faa"
outname = arguments[3]
#outname = "/home/nmoreyra/Data/cbuzz_genomes/09_detox-genes/drosophila_datasets/gff_faa_files/Ccap.1prot.faa"

#species_list = ['Ccap','Dana','Dari','Dato','Dbrb','Dbuz','Dele','Dere','Dficu','Dgri', 'Dgua','Dhyd','DkoeA','DkoeB','Dmau','Dmela','Dmex','Dmoj','Dnav','Dnova','Dobs','Dper','Dpse','Drho','Dsec','Dsim','Dsub','Dsuz','Dtaka','Dvir','Dwil','Dyak','Gfus','Mdom']
genes2nprots, gene2prot = get_genes_multiples_prots(gff)
prot_seqs = dict_pep(prots, biopython=True)

with open(outname, 'wt') as OUT:
	for gene,n in genes2nprots.items():
		if n > 1:
			prot_list = gene2prot[gene]
			longest_prot_id = find_longest_prot(prot_list, prot_seqs)
			OUT.write(">" + longest_prot_id + "\n" + prot_seqs[longest_prot_id] + "\n")
		else:
			protid = gene2prot[gene][0]
			OUT.write(">" + protid + "\n" + prot_seqs[protid] + "\n")
# end of main program













