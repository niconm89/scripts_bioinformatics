#!/usr/bin/env python3

"""
Script: fasta_basic_stats.py
Author: Nicolás Nahuel Moreyra (niconm89@gmail.com)
GitHub: @niconm89
Date: Mon Feb 12 21:46:01 2024

Description:
    This script calculates statistics for a FASTA file, including the species name, number of sequences,
    total length of sequences, and mean length of sequences.

Usage:
    python fasta_basic_stats.py

Requirements:
    - Biopython library (https://biopython.org/)
"""

#%% Imports
from Bio import SeqIO
import os
import gzip

#%% Functions
def fasta_stats(fasta_file):
    """
    Calculate statistics for a FASTA file.

    Args:
        fasta_file (str): Path to the input FASTA file.

    Returns:
        tuple: A tuple containing the species name, number of sequences, total length of sequences,
               and mean length of sequences.
    """
    # Extract the species name from the file name
    species_name = os.path.basename(fasta_file).split('.')[0]

    # Initialize variables for sequence count and total length
    num_sequences = 0
    total_length = 0

    # Open the FASTA file, compressed or uncompressed
    if fasta_file.endswith('.gz'):
        handle = gzip.open(fasta_file, "rt")
    else:
        handle = open(fasta_file, "r")

    # Iterate over sequences in the FASTA file
    for record in SeqIO.parse(handle, "fasta"):
        num_sequences += 1
        total_length += len(record.seq)

    # Calculate the mean length of sequences
    mean_length = total_length / num_sequences if num_sequences > 0 else 0

    return species_name, num_sequences, total_length, mean_length

if __name__ == "__main__":
    # Prompt user for the input FASTA file
    fasta_file = input("Enter the path to the input FASTA file: ")

    # Check if the FASTA file exists
    if not os.path.exists(fasta_file):
        print("Error: File '{}' not found.".format(fasta_file))
        sys.exit(1)

    # Calculate statistics for the FASTA file
    species_name, num_sequences, total_length, mean_length = fasta_stats(fasta_file)

    # Prompt user for the output file name
    output_file = input("Enter the name for the output file to save the results: ")

    # Write the results to the output file
    with open(output_file, "w") as f:
        f.write("Species: {}\n".format(species_name))
        f.write("Number of sequences: {}\n".format(num_sequences))
        f.write("Total length of sequences: {}\n".format(total_length))
        f.write("Mean length of sequences: {}\n".format(mean_length))

    print("Statistics saved to '{}'.".format(output_file))
