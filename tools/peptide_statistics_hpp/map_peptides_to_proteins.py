#!/usr/bin/python


import sys
import getopt
import os
import json
import ming_protein_library
import ming_fileio_library
import re
from collections import defaultdict


def usage():
    print("<input fasta> <input peptides> <output results file>")

def strip_sequence(sequence):
    if sequence[1] == "." and sequence[-2] == ".":
        sequence = sequence[2:-2]

    p = re.compile('\W|\d')
    sequence = p.sub("", sequence)
    return sequence

def load_peptide_list(input_peptide_filename):
    row_count, table_data = ming_fileio_library.parse_table_with_headers(input_peptide_filename)
    all_peptides = []
    if "Peptide" in table_data:
        all_peptides = table_data["Peptide"]
    elif "peptide" in table_data:
        all_peptides = table_data["peptide"]
    elif "Annotation" in table_data:
        all_peptides = table_data["Annotation"]
    elif "Sequence" in table_data:
        all_peptides = table_data["Sequence"]
    elif "variant_sequence" in table_data:
        all_peptides = table_data["variant_sequence"]
    elif "sequence" in table_data:
        all_peptides = table_data["sequence"]


    return_peptides = []

    for peptide in all_peptides:
        return_peptides.append(peptide)


    return_peptides = list(set(return_peptides))
    return_peptides.sort()

    return return_peptides

def main():
    input_fasta_filename = sys.argv[1]
    input_peptide_list_filename = sys.argv[2]
    output_coverage_filename = sys.argv[3]

    proteome = ming_protein_library.parse_fasta_proteome_file(input_fasta_filename)
    #IL substitution for proteome
    for protein in proteome.protein_list:
        protein.sequence = protein.sequence.replace("I", "L")

    peptide_list = load_peptide_list(input_peptide_list_filename)
    protein_list = proteome.protein_list

    #IL peptides
    il_peptide_list = list(set([strip_sequence(x).replace("I", "L") for x in peptide_list]))

    peptide_to_protein_map = proteome.get_peptides_mapped_to_proteins_efficient_w_coordinates(il_peptide_list)

    output_dict = defaultdict(list)

    for peptide in peptide_list:
        proteins_covered_w_coords = peptide_to_protein_map[strip_sequence(peptide).replace("I", "L")]
        proteins_covered = list(set([p.split("(")[0] for p in proteins_covered_w_coords]))
        output_dict["peptide"].append(peptide)
        output_dict["il_peptide"].append(strip_sequence(peptide).replace("I", "L"))
        output_dict["proteins"].append(";".join(proteins_covered))
        output_dict["proteins_w_coords"].append(";".join(proteins_covered_w_coords))
        output_dict["num_proteins"].append(len(proteins_covered))

    print("Final output")
    ming_fileio_library.write_dictionary_table_data(output_dict, output_coverage_filename)


if __name__ == "__main__":
    main()
