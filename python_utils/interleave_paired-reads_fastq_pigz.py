#!/usr/bin/env python3

"""
Script: interleave_paired-reads_fastq_pigz.py
Author: Nicolás Nahuel Moreyra (niconm89@gmail.com)
GitHub: @niconm89
Date: Mon Apr 22 14:11:46 2024

Description:
    
Usage:
    python interleave_paired-reads_fastq_pigz.py

Requirements:
    - subprocess, pigz (on system and available in $PATH)
"""

#%% Imports
import subprocess

#%% Functions
def pigz_compress(input_file, output_file):
    """
    Compress a file using pigz.
    """
    subprocess.run(["pigz", "-9", "-c", input_file, "-o", output_file])


def pigz_decompress(input_file, output_file):
    """
    Decompress a file using pigz.
    """
    subprocess.run(["pigz", "-d", "-c", input_file, "-o", output_file])
#end
def interleaver(forward_lines, reverse_lines):
    interleaved_lines = []
    for f_line, r_line in zip(forward_lines, reverse_lines):
        interleaved_lines.extend([f_line, r_line])
    return interleaved_lines
#end

def interleave_reads(forward_file, reverse_file, interleaved_output_file):
    """
    Interleave paired reads from forward and reverse FASTQ files into a single interleaved FASTQ file.
    """
    with subprocess.Popen(["pigz", "-dc", forward_file], stdout=subprocess.PIPE) as forward, \
            subprocess.Popen(["pigz", "-dc", reverse_file], stdout=subprocess.PIPE) as reverse, \
            subprocess.Popen(["pigz", "-c"], stdin=subprocess.PIPE, stdout=open(interleaved_output_file, 'wb')) as interleaved: #open(interleaved_file, 'wb') as interleaved:

        # Define chunk size for buffering reads
        chunk_size = 2000000  # Adjust as needed 

        while True:
            forward_chunk = [forward.stdout.readline() for _ in range(4 * chunk_size)]
            reverse_chunk = [reverse.stdout.readline() for _ in range(4 * chunk_size)]
            if not forward_chunk[0] or not reverse_chunk[0]:
                break

            interleaved_chunk = interleaver(forward_chunk, reverse_chunk)
            interleaved.stdin.writelines(interleaved_chunk)
#end

if __name__ == '__main__':
    forward_reads_file = input("Enter the path to the forward reads FASTQ file: ")
    reverse_reads_file = input("Enter the path to the reverse reads FASTQ file: ")
    interleaved_file = input("Enter the name for the output interleaved file: ")

    # Call the function to interleave reads
    if interleaved_file.endswith("gz"):
        interleaved_output_file = interleaved_file
    else:
        interleaved_output_file = interleaved_file + ".gz"
    interleave_reads(forward_reads_file, reverse_reads_file, interleaved_output_file)

    print(f"Interleaved reads have been written to {interleaved_file}.gz")
