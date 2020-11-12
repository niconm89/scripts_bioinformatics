#!/usr/bin/python
#NUcleotide Diversity Estimation : NUDE

import sys
import getopt
from gendivest import gde

def usage():
	print """
nude.py : reads a clustalw alignment file and builds a file with the nucleotide diversity values (windows) between the sequences

nude.py [-h] [-al <length>] <filename>

-h		print this message
-al <alignment>	clustalw alignment file


<filename>	the file has to be in clustalw2 format

	"""


o,a = getopt.getopt(sys.argv[1:],'l:h') #o es una lista de argumentos opcionales y a es una lista de elementos requeridos


#dna=raw_input("Enter a DNA sequence, please: ") #input() para python3

filename=sys.argv[1]

readingframe = int(sys.argv[2])

id_sequence = sys.argv[3]

try: 
	f = open(filename)
except IOError:
	print("File %s does not exist!!" % filename)
	
##Largo de las secuencias
#from functions_exam import lengths
#lista_longitudes = lengths(records)

##Secuencia más corta
#from functions_exam import shortest_seq
#minimo,rep2,identificadores2 = shortest_seq(records)

"""

here I can write whatever I want, this is a comments space!

"""

#print ("The gc content is %8.3f%%" % gc_percent)


