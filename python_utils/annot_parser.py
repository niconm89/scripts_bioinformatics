#!/usr/bin/env python3

"""
Script: annot_parser.py
Author: Nicolás Nahuel Moreyra (niconm89@gmail.com)
GitHub: @niconm89
Date: Wed Oct 04 15:03:10 2023

Description:
    

Usage:
    python annot_parser.py -h

Requirements:
    - argparse
"""

#%% Imports
import argparse

#%% Functions
def parse_gff(file_path, feature, pattern, attributes_to_print):
    """
    Parsea un archivo GFF en busca de patrones y atributos específicos.

    Args:
        file_path (str): Ruta al archivo GFF.
        patterns (list of str): Lista de patrones a buscar en la descripción.
        attributes_to_print (list of str): Lista de atributos a imprimir.

    Returns:
        None
    """
    with open(file_path, 'r') as f:
        for line in f:
            if not line.startswith('#'):  # Ignora las líneas de comentarios
                fields = line.strip().split('\t')
                if len(fields) >= 9 and fields[2] == feature:  # Verifica si es una línea de gene
                    attributes_gene = dict(item.split('=') for item in fields[8].split(';'))
                    description_gene = switch_case(feature, attributes_gene)
                    # Busca el patron en la descripción del gen
                    if pattern.lower() in description_gene.lower():
                        #print(attributes_gene)
                        # Imprime los atributos especificados
                        for attrr in attributes_to_print:
                            if attrr in attributes_gene:
                                print(f"{attrr.capitalize()}: {attributes_gene[attrr]}")

def switch_case(feature, attributes_gene):
    if feature == "gene":
        description = attributes_gene.get('description', 'N/A')  # Obtiene la descripción o 'N/A' si no existe
    elif feature == "mRNA":
        description = attributes_gene.get('product', 'N/A')  # Obtiene la descripción o 'N/A' si no existe
    elif feature == "exon":
        description = attributes_gene.get('product', 'N/A')  # Obtiene la descripción o 'N/A' si no existe
    elif feature == "CDS":
        description = attributes_gene.get('product', 'N/A')  # Obtiene la descripción o 'N/A' si no existe
    else:
        print("Invalid case")
        exit(1)
    return description


def main():
    """
    Función principal para manejar los argumentos de línea de comandos y llamar a parse_gff.
    """
    parser = argparse.ArgumentParser(description="Parse GFF file")
    parser.add_argument("file_path", help="Path to the GFF file")
    parser.add_argument("feature", help="Patterns to search for")
    parser.add_argument("pattern", help="Patterns to search for")
    parser.add_argument("--attributes", nargs='+', default=["ID","description"], help="Attributes to print. By defult just de gene_id and its description are shown.")
    args = parser.parse_args()

    # Llama a la función parse_gff con los argumentos proporcionados
    parse_gff(args.file_path, args.feature, args.pattern, args.attributes)

if __name__ == "__main__":
    main()
