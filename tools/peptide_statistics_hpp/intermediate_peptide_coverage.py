#!/usr/bin/python


import sys
import getopt
import os
import json
import ming_protein_library
import ming_fileio_library
import re
import cPickle as pickle
from collections import defaultdict

def usage():
    print("<input fasta> <input peptides> <input parameters file> <output results folder>")

# def strip_sequence(sequence):
#     if sequence[1] == "." and sequence[-1] == ".":
#         sequence = sequence[1:-1]
#
#     p = re.compile('\W|\d')
#     sequence = p.sub("", sequence)
#     return sequence

def strip_sequence(sequence):
    if sequence[1] == "." and sequence[-2] == ".":
        sequence = sequence[2:-2]
    # whitelist =  ["57.021464","15.994915","14.015650","17.026549","43.005814","42.010565","0.984016"]
    # whitelist += ["57.021",   "15.995",   "14.016",   "17.027",   "43.006",   "42.011",   "0.984"   ,"229.163"]

    # for mod in whitelist:
    #     sequence = sequence.replace(mod,"")
    # if len([s for s in sequence if s.isdigit()]) == 0:
    p = re.compile('\W|\d')
    sequence = p.sub("", sequence)
    # else:
    #     sequence = None
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
        stripped_peptide = strip_sequence(peptide)
        return_peptides.append(stripped_peptide)

    return_peptides = list(set(return_peptides))

    return return_peptides

def main():
    input_fasta_filename = sys.argv[1]
    input_peptide_list_filename = sys.argv[2]
    input_parameters_filename = sys.argv[3]
    input_nextprot_pe = sys.argv[4]
    output_results_folder = sys.argv[5]

    partition_total = 1
    partition_of_node = 0

    print("LOADING " + input_parameters_filename)

    params_map = json.loads(open(input_parameters_filename).read())
    partition_total = params_map["total_paritions"]
    partition_of_node = params_map["node_partition"]

    output_coverage_filename = os.path.join(output_results_folder, str(partition_of_node) + ".json")

    proteome = ming_protein_library.parse_fasta_proteome_file(input_fasta_filename,input_nextprot_pe)
    #IL substitution for proteome
    for protein in proteome.protein_list:
        protein.sequence = protein.sequence.replace("I", "L")

    peptide_list = load_peptide_list(input_peptide_list_filename)
    protein_list = proteome.protein_list

    #Lets alternate what we get
    peptide_list = list(set([strip_sequence(x) for x in peptide_list[partition_of_node::partition_total] if strip_sequence(x) != None]))

    #IL peptides
    il_peptide_list = list(set([x.replace("I", "L") for x in peptide_list]))

    protein_to_unique_peptide_map = defaultdict(list)
    peptide_to_protein_map = proteome.get_peptides_mapped_to_proteins_efficient_w_coordinates(il_peptide_list)

    for peptide in peptide_list:
        if peptide:
            proteins_covered = peptide_to_protein_map[peptide.replace("I", "L")]
            print(proteins_covered)




if __name__ == "__main__":
    main()
