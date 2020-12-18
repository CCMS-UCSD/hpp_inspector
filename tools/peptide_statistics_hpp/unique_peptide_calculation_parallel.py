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
    whitelist =  ["57.021464","15.994915","14.015650","17.026549","43.005814","42.010565","0.984016"]
    whitelist += ["57.021",   "15.995",   "14.016",   "17.027",   "43.006",   "42.011",   "0.984"   ,"229.163"]

    for mod in whitelist:
        sequence = sequence.replace(mod,"")
    if len([s for s in sequence if s.isdigit()]) == 0:
        p = re.compile('\W|\d')
        sequence = p.sub("", sequence)
    else:
        sequence = None
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
        if stripped_peptide:
            stripped_peptide = stripped_peptide.replace("I", "L")
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
    proteome.create_decoy_proteins()
    
    protein_to_gene = {}


    #IL substitution for proteome
    for protein in proteome.protein_list:
        protein.sequence = protein.sequence.replace("I", "L")
        protein_to_gene[protein.protein] = protein.gene_name

    peptide_list = load_peptide_list(input_peptide_list_filename)
    protein_list = proteome.protein_list

    #Lets alternate what we get
    peptide_list = peptide_list[partition_of_node::partition_total]

    protein_to_gene_unique_peptide_map_no_mismatch = defaultdict(list)
    protein_to_gene_unique_peptide_map = defaultdict(list)
    protein_to_protein_unique_peptide_map_no_mismatch = defaultdict(list)
    protein_to_protein_unique_peptide_map = defaultdict(list)
    peptide_list_gte9 = [p for p in peptide_list if len(p) >= 9]
    peptide_to_protein_map_w_mismatch = proteome.get_peptides_mapped_to_proteins_mismatched(peptide_list_gte9, mismatch_max=1)
    peptide_to_protein_map = proteome.get_peptides_mapped_to_proteins_efficient_w_coordinates(peptide_list)

    for peptide in peptide_list:
        shared_gene = True
        min_genes_covered = None
        for p in peptide_to_protein_map_w_mismatch[peptide]:
            genes_covered = protein_to_gene[p]
            if min_genes_covered is None:
                min_genes_covered = genes_covered
            else:
                shared_with_protein = False
                if len(min_genes_covered) < len(genes_covered) and genes_covered[len(min_genes_covered)] is '-':
                    if min_genes_covered in genes_covered:
                        shared_with_protein = True
                elif len(min_genes_covered) > len(genes_covered) and min_genes_covered[len(genes_covered)] is '-':
                    if genes_covered in min_genes_covered:
                        min_genes_covered = genes_covered
                        shared_with_protein = True
                elif len(min_genes_covered) == len(genes_covered):
                    if genes_covered in min_genes_covered:
                        shared_with_protein = True
                shared_gene = shared_gene and shared_with_protein
        if shared_gene and len(peptide) >= 9 and min_genes_covered:
            protein_pos_covered = peptide_to_protein_map[peptide]
            for protein_pos in protein_pos_covered:
                protein,both_pos = protein_pos.split(" ")
                pos = both_pos.replace("(","").replace(")","").split('-')
                start = int(pos[0])
                end = int(pos[1])
                protein_to_gene_unique_peptide_map[protein].append((peptide, (start,end)))

        shared_gene = True
        min_genes_covered = None
        for protein_pos in peptide_to_protein_map[peptide]:
            p,both_pos = protein_pos.split(" ")
            genes_covered = protein_to_gene[p]
            if min_genes_covered is None:
                min_genes_covered = genes_covered
            else:
                shared_with_protein = False
                if len(min_genes_covered) < len(genes_covered) and genes_covered[len(min_genes_covered)] is '-':
                    if min_genes_covered in genes_covered:
                        shared_with_protein = True
                elif len(min_genes_covered) > len(genes_covered) and min_genes_covered[len(genes_covered)] is '-':
                    if genes_covered in min_genes_covered:
                        min_genes_covered = genes_covered
                        shared_with_protein = True
                elif len(min_genes_covered) == len(genes_covered):
                    if genes_covered in min_genes_covered:
                        shared_with_protein = True
                shared_gene = shared_gene and shared_with_protein
        if shared_gene and len(peptide) >= 9 and min_genes_covered:
            protein_pos_covered = peptide_to_protein_map[peptide]
            for protein_pos in protein_pos_covered:
                protein,both_pos = protein_pos.split(" ")
                pos = both_pos.replace("(","").replace(")","").split('-')
                start = int(pos[0])
                end = int(pos[1])
                protein_to_gene_unique_peptide_map_no_mismatch[protein].append((peptide, (start,end)))

        protein_pos_covered = peptide_to_protein_map[peptide]
        if len(set(protein_pos.split(" ")[0] for protein_pos in protein_pos_covered)) == 1:
            for protein_pos in protein_pos_covered:
                protein,both_pos = protein_pos.split(" ")
                pos = both_pos.replace("(","").replace(")","").split('-')
                start = int(pos[0])
                end = int(pos[1])
                protein_to_protein_unique_peptide_map_no_mismatch[protein].append((peptide, (start,end)))

        protein_pos_covered = peptide_to_protein_map[peptide]
        if len(set(peptide_to_protein_map_w_mismatch[peptide])) == 1:
            for protein_pos in protein_pos_covered:
                protein,both_pos = protein_pos.split(" ")
                pos = both_pos.replace("(","").replace(")","").split('-')
                start = int(pos[0])
                end = int(pos[1])
                protein_to_protein_unique_peptide_map[protein].append((peptide, (start,end)))


    output_dict = {}
    output_dict["protein_to_protein_unique_peptide_no_mismatch"] = protein_to_protein_unique_peptide_map_no_mismatch
    output_dict["protein_to_protein_unique_peptide"] = protein_to_protein_unique_peptide_map
    output_dict["protein_to_gene_unique_peptide_no_mismatch"] = protein_to_gene_unique_peptide_map_no_mismatch
    output_dict["protein_to_gene_unique_peptide"] = protein_to_gene_unique_peptide_map
    output_dict["peptide_to_protein_list"] = peptide_to_protein_map

    open(output_coverage_filename, "w").write(json.dumps(output_dict))


if __name__ == "__main__":
    main()
