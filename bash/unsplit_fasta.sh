#!/usr/bin/env bash

"""
Script: unsplit_fasta.sh
#Author: Nicolás Nahuel Moreyra (niconm89@gmail.com)
#GitHub: @niconm89
#Date: Mon Apr 22 10:36:12 2024

Description:
    This script convert a FASTA file containing multiple lines per sequence into a one line per sequence (plus de sequence ID line) FASTA file.
    
Usage:
    bash unsplit_fasta.sh fasta_file_name

Requirements:
    - awk
"""

awk '/^>/ {printf("%s%s\n",(NR>1)?"\n":"" ,$0);next; } { printf("%s",$0);} END {printf("\n");}' $1
