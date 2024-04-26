#!/usr/bin/env python3


"""
Script: filter_fasta_by_length.py
Author: Nicolás Nahuel Moreyra (niconm89@gmail.com)
GitHub: @niconm89
Date: Tue Nov 07 12:11:56 2023

Description:
    
Usage:
    python filter_fasta_by_length.py

Requirements:
    - Biopython
"""

#%% Import
from Bio import SeqIO

#%% Functions
def filter_fasta(input_file, output_file, min_length, max_length):
    """
    Filtra un archivo FASTA y guarda las secuencias cuya longitud está dentro del rango especificado.

    Parámetros:
        input_file (str): Nombre del archivo FASTA de entrada.
        output_file (str): Nombre del archivo FASTA de salida filtrado.
        min_length (int): Longitud mínima deseada para las secuencias.
        max_length (int): Longitud máxima deseada para las secuencias.
    """
    with open(output_file, 'w') as out_fasta:
        for record in SeqIO.parse(input_file, 'fasta'):
            seq_length = len(record.seq)
            if min_length <= seq_length <= max_length:
                SeqIO.write(record, out_fasta, 'fasta')

def main():
    print("Bienvenido al filtro de archivos FASTA")
    input_file = input("Ingrese el nombre del archivo FASTA de entrada: ")
    output_file = input("Ingrese el nombre del archivo FASTA de salida filtrado: ")
    min_length = int(input("Ingrese la longitud mínima deseada para las secuencias: "))
    max_length = int(input("Ingrese la longitud máxima deseada para las secuencias: "))

    filter_fasta(input_file, output_file, min_length, max_length)
    print(f"Se han filtrado las secuencias y guardado en '{output_file}'.")

if __name__ == "__main__":
    main()
