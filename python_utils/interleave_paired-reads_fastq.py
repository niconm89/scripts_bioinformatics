#!/usr/bin/env python3

"""
Script: interleave_paired-reads_fastq.py
Author: Nicolás Nahuel Moreyra (niconm89@gmail.com)
GitHub: @niconm89
Date: Mon Apr 22 14:11:46 2024

Description:
    
Usage:
    python interleave_paired-reads_fastq.py

Requirements:
    - gzip
"""

#%% Imports
import gzip

#%% Functions
def interleave_reads(forward_file, reverse_file, interleaved_file):
    """
    Interleave paired reads from forward and reverse FASTQ files into a single interleaved FASTQ file.
    """
    with gzip.open(forward_file, 'rt') as forward, \
            gzip.open(reverse_file, 'rt') as reverse, \
            gzip.open(interleaved_file, 'wt') as interleaved:
        while True:
            # Read four lines at a time from both forward and reverse files
            forward_lines = [forward.readline().strip() for _ in range(4)]
            reverse_lines = [reverse.readline().strip() for _ in range(4)]

            # If any of the lines are empty, we've reached the end of the file
            if not forward_lines[0] or not reverse_lines[0]:
                break

            # Write interleaved reads to the output file
            interleaved.write('\n'.join(forward_lines) + '\n')
            interleaved.write('\n'.join(reverse_lines) + '\n')


if __name__ == '__main__':
    # Prompt user for the input file names
    forward_reads_file = input("Enter the path to the forward reads FASTQ file: ")
    reverse_reads_file = input("Enter the path to the reverse reads FASTQ file: ")
    interleaved_file = input("Enter the name for the output interleaved file: ")

    # Call the function to interleave reads
    interleave_reads(forward_reads_file, reverse_reads_file, interleaved_file)

    print(f"Interleaved reads have been written to {interleaved_file}")
