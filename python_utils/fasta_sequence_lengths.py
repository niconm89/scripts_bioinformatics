#!/usr/bin/env python3

"""
Script: fasta_sequence_lengths.py
Author: Nicolás Nahuel Moreyra (niconm89@gmail.com)
GitHub: @niconm89
Date: Fri Apr 19 12:31:27 2024

Description:
    
Usage:
    python fasta_sequence_lengths.py

Requirements:
    - 
"""
#%% Imports

#%% Functions
def read_fasta(file_path):
    """
    Read a FASTA file and yield each sequence identifier and its corresponding sequence.
    """
    with open(file_path, 'r') as fasta_file:
        sequence_id = None
        sequence = ''
        for line in fasta_file:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if sequence_id:
                    yield sequence_id, sequence
                sequence_id = line[1:]
                sequence = ''
            else:
                sequence += line
        if sequence_id:
            yield sequence_id, sequence


def write_sequence_lengths(fasta_file_path, output_file_path):
    """
    Read a FASTA file, compute the length of each sequence, and write the lengths into a text file.
    """
    with open(output_file_path, 'w') as output_file:
        for seq_id, seq in read_fasta(fasta_file_path):
            length = len(seq)
            output_file.write(f'{seq_id}\t{length}\n')


def main():
    print("Welcome to the FASTA sequence length calculator!")
    input_fasta = input("Please enter the path to the input FASTA file: ")
    output_file = input("Please enter the name for the output file: ")

    # Call the function to write sequence lengths to the output file
    write_sequence_lengths(input_fasta, output_file)

    print(f"Sequence lengths have been written to {output_file}")


if __name__ == '__main__':
    main()
