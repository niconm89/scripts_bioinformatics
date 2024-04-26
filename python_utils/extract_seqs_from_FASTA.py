#!/usr/bin/env python3

"""
Script: extract_seqs_from_FASTA.py
Author: Nicolás Nahuel Moreyra (niconm89@gmail.com)

Date: Tue Mar 17 09:11:56 2021

Description:
    This script extracts sequences from a FASTA file based on different extraction modes:
    - head: extracts the first n sequences.
    - tail: extracts the last n sequences.
    - random: randomly selects n sequences.
    - byIDs: extracts sequences based on a list of sequence IDs provided in a separate file.

Usage:
    python extract_seqs_from_FASTA.py -f <input_file> -o <output_file> -m <mode> [-n <nseqs>] [-d <IDs_file>]

Arguments:
    -f, --infile: Full path to the input FASTA file to extract sequences from. It can be gzipped.
    -o, --outfile: Full path to the output file where the extracted sequences will be saved. Output file will not be compressed.
    -m, --mode: Extraction mode. Choose one of the following options: 'head', 'tail', 'random', 'byIDs'.
    -n, --nseqs: Number of sequences to recover from the input. Required for 'head', 'tail', and 'random' modes.
    -d, --IDsfile: Full path to the file containing sequence IDs to extract. Required for 'byIDs' mode. Sequence descriptions must be excluded.

Requirements:
    - Biopython library (https://biopython.org/)
"""

#%% Imports
import argparse
import os
from Bio import SeqIO
import random # comment this line if you are not using the random extraction mode and don't want to install the library
import gzip

#%% Function definitions
def check_paths(inputfile,outputfile):
    """
    Check if the input file and output folder exist.

    Args:
        inputfile (str): Path to the input FASTA file.
        outputfile (str): Path to the output file.

    Raises:
        FileNotFoundError: If the input file or output folder does not exist.
    """
    if not os.path.isfile(inputfile):
        raise FileNotFoundError("Error: Input file does not exist! " + inputfile)
    
    if '/' in outputfile:
        elements_path = outputfile.split('/')
        outputfolder = '/'.join(elements_path[:-1])
        if not os.path.isdir(outputfolder):
            raise FileNotFoundError("Error: Output folder does not exist! " + outputfolder)
#end
def check_mode(arguments):
    """
    Check if the chosen extraction mode is valid and all required parameters are provided.

    Args:
        arguments (argparse.Namespace): Parsed arguments from command line.

    Raises:
        ValueError: If the extraction mode is invalid or required parameters are missing.
    """
    if arguments.mode == 'byIDs':
        if arguments.IDsfile:
            if not os.path.isfile(arguments.IDsfile):
                raise FileNotFoundError("Error: IDs file does not exist! " + arguments.IDsfile)
        else:
            raise ValueError("Extraction mode was set to byIDs but no IDs file was input. Aborting...")
    else:
        if not arguments.nseqs:
            raise ValueError("Unless you are using byIDs extraction mode, you must define the number of sequences you want to recover (-n, --nseqs). Aborting...")
#end
def read_IDs(IDsfile):
    """
    Read sequence IDs from a file and return them as a list.

    Args:
        IDsfile (str): Path to the file containing sequence IDs.

    Returns:
        list: List of sequence IDs.
    """
    IDs_list = []
    with open(IDsfile,'rt') as inputIDs:
        lines = inputIDs.readlines()
        for line in lines:
            id_seq = line.split(" ")[0]
            IDs_list.append(id_seq.rstrip())
    return IDs_list
#
def check_gzipped(inputfile):
    """
    Check if a file is gzipped.

    Args:
        inputfile (str): Path to the input file.

    Returns:
        bool: True if the file is gzipped, False otherwise.
    """
    # The first two bytes of a gzip file are: 1f8b
    gzipped = False
    GZIP_MAGIC_NUMBER = b'\x1f\x8b'
    with open(inputfile,'rb') as INPUT:
        twobytes = INPUT.read(2)
        if twobytes == GZIP_MAGIC_NUMBER:
            gzipped = True
    return gzipped
#end
def read_fasta(inputfile):
    """
    Read sequences from a FASTA file.

    Args:
        inputfile (str): Path to the input FASTA file.

    Returns:
        dict: Dictionary containing sequence IDs as keys and sequences as values.
    """
    sequences = {}
    if check_gzipped(inputfile):
        with gzip.open(inputfile, "rt") as FASTA:
            count = 0
            for record in SeqIO.parse(FASTA, "fasta"):
                count +=1
                seqname = str(record.id)
                seqprot = str(record.seq)
                if seqname in sequences:
                    print("There are repeated sequence IDs in the input file. Aborting...")
                    exit(0)
                sequences[seqname] = seqprot
    else:
        with open(inputfile, "rt") as FASTA:
            for record in SeqIO.parse(FASTA, "fasta"):
                seqname = str(record.id)
                seqprot = str(record.seq)
                if seqname in sequences:
                    print("There are repeated sequence IDs in the input file. Aborting...")
                    exit(0)
                sequences[seqname] = seqprot
    return sequences
#end
def extract_nseqs(inputfile, outputfile, nseqs, mode, IDsfile):
    """
    Extract sequences from the input FASTA file based on the chosen extraction mode.

    Args:
        inputfile (str): Path to the input FASTA file.
        outputfile (str): Path to the output file to save the extracted sequences.
        nseqs (int): Number of sequences to extract (required for 'head', 'tail', and 'random' modes).
        mode (str): Extraction mode ('head', 'tail', 'random', 'byIDs').
        IDsfile (str): Path to the file containing sequence IDs (required for 'byIDs' mode).
    """
    sequences = read_fasta(inputfile)
    subseqs = {}
    if mode == 'head' or mode == 'tail':
        subseqs = extract_head_tail(sequences, nseqs, mode)
    else:
        if mode == 'random':
            subseqs = extract_random(sequences, nseqs)
        else:
            IDs_list = read_IDs(IDsfile)
            subseqs = extract_byIDs(sequences, IDs_list)
    save_subseqs(subseqs, outputfile)
#end
def save_subseqs(subseqs, outputfile):
    """
    Save extracted sequences to the output file.

    Args:
        subseqs (dict): Dictionary containing sequence IDs as keys and sequences as values.
        outputfile (str): Path to the output file.
    """
    with open(outputfile, 'wt') as OUT:
        for k,v in subseqs.items():
            OUT.write(">" + str(k) + "\n" + str(v) + "\n")
#end
def extract_head_tail(sequences, nseqs, mode):
    """
    Extract sequences from the beginning or end of the input sequence dictionary.

    Args:
        sequences (dict): Dictionary containing sequence IDs as keys and sequences as values.
        nseqs (int): Number of sequences to extract.
        mode (str): Extraction mode ('head' or 'tail').

    Returns:
        dict: Dictionary containing extracted sequences.
    """
    subsequences = {}
    if mode == 'head':
        for ID in list(sequences)[0:nseqs]:
            subsequences[ID] = sequences[ID]
    else:
        pos = len(sequences) - nseqs
        for ID in list(sequences)[pos:]:
            subsequences[ID] = sequences[ID]
    return subsequences
#end
def extract_random(sequences, nseqs):
    """
    Extract random sequences from the input sequence dictionary.

    Args:
        sequences (dict): Dictionary containing sequence IDs as keys and sequences as values.
        nseqs (int): Number of sequences to extract.

    Returns:
        dict: Dictionary containing extracted sequences.
    """
    subsequences = {}
    lista_IDs = list(sequences.keys())
    for n in range(0, nseqs):
        IDseq = random.choice(lista_IDs)
        lista_IDs.remove(IDseq)
        subsequences[IDseq] = sequences[IDseq]
    return subsequences
#end
def extract_byIDs(sequences, IDs_list):
    """
    Extract sequences based on a list of sequence IDs.

    Args:
        sequences (dict): Dictionary containing sequence IDs as keys and sequences as values.
        IDs_list (list): List of sequence IDs to extract.

    Returns:
        dict: Dictionary containing extracted sequences.
    """
    subsequences = {}
    for ID in IDs_list:
        if ID in sequences:
            subsequences[ID] = sequences[ID]
        else:
            print("Sequence ID " + str(ID) + " not found. ID avoided...")
    return subsequences
#end

#%%
def main():
    """
    Main function to parse command line arguments and execute sequence extraction.
    """
    parser = argparse.ArgumentParser(
        description='extract_seqs_from_FASTA.py recover sequences from a FASTA file.',
        epilog="""End of the help""")

    parser.add_argument('-f', '--infile', type=str, required=True, help='Full path to the FASTA file to convert. It can be gzipped.')
    parser.add_argument('-o', '--outfile', type=str, required=True, help='Full path where the new FASTA file will be placed. Output file will not be compressed.')
    parser.add_argument('-m', '--mode', type=str, required=True, choices=['head', 'tail', 'random', 'byIDs'], help='Extraction mode. head: first -n sequences; tail: last -n sequences; random: randomly selection of -n sequences; ID_file: file containing a sequence ID per line.')
    parser.add_argument('-n', '--nseqs', type=int, required=False, help='Numbers of sequences to recover from the input.')
    parser.add_argument('-d', '--IDsfile', type=str, required=False, help='Full path to the file with IDs to extract. Sequence descriptions must be excluded.')

    args=parser.parse_args()

    check_paths(args.infile, args.outfile) # check if input file and output folder exist
    check_mode(args)

    extract_nseqs(args.infile, args.outfile, args.nseqs, args.mode, args.IDsfile)

#%% Main program
if __name__ == '__main__':
    main()
