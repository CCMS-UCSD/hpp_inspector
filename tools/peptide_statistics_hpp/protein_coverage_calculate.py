#!/usr/bin/python


import sys
import getopt
import os
import json
import ming_protein_library
import uniprot_domain_statistics
import ming_fileio_library
import re
import ming_numerical_utilities

if (sys.version_info > (3, 0)):
    # Python 3 code in this block
    import _pickle as pickle
else:
    # Python 2 code in this block
    import cPickle as pickle

def usage():
    print("<input fasta> <input peptides> <input parameters file> <input gff pkl file> <output results folder>")

# def strip_sequence(sequence):
#     if sequence[2] == "." and sequence[-2] == ".":
#         sequence = sequence[2:-2]
#     # whitelist =  ["57.021464","15.994915","14.015650","17.026549","43.005814","42.010565","0.984016"]
#     # whitelist += ["57.021",   "15.995",   "14.016",   "17.027",   "43.006",   "42.011",   "0.984"   ,"229.163"]
#
#     # for mod in whitelist:
#     #     sequence = sequence.replace(mod,"")
#     # if len([s for s in sequence if s.isdigit()]) == 0:
#     p = re.compile('\W|\d')
#     sequence = p.sub("", sequence)
#     # else:
#     #     sequence = None
#     return sequence

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
    input_uniprot_gff_pkl_filename = sys.argv[4]
    input_nextprot_pe = sys.argv[5]
    output_results_folder = sys.argv[6]

    output_protein_coverage_folder = sys.argv[7]

    partition_total = 1
    partition_of_node = 0

    # print("LOADING " + input_parameters_filename)

    params_map = json.loads(open(input_parameters_filename).read())
    partition_total = params_map["total_paritions"]
    partition_of_node = params_map["node_partition"]

    output_coverage_filename = os.path.join(output_results_folder, str(partition_of_node) + ".json")

    proteome = ming_protein_library.parse_fasta_proteome_file(input_fasta_filename,input_nextprot_pe)
    proteome.create_decoy_proteins()

    protein_list = proteome.protein_list
    protein_count = len(protein_list)
    #Lets alternate what we get
    protein_list = protein_list[partition_of_node::partition_total]
    proteome.protein_list = protein_list

    #IL substitution for proteome
    for protein in proteome.protein_list:
        protein.sequence = protein.sequence.replace("I", "L")

    #protein_to_annotations_map = json.load(open(input_uniprot_gff_pkl_filename, "r"))
    print("Loading PKL of Uniprot annotations")
    protein_to_annotations_map = pickle.load(open(input_uniprot_gff_pkl_filename, "rb"))

    # peptide_list = load_peptide_list(input_peptide_list_filename)
    peptide_list = list(set([strip_sequence(x) for x in load_peptide_list(input_peptide_list_filename) if strip_sequence(x) != None]))

    protein_coverage_map = {}

    print("Mapping peptides to proteins")
    proteins_to_peptide_map = proteome.get_proteins_to_peptides_efficient(peptide_list)

    debug_coverage_dict = {}

    progress_count = 0
    for protein in protein_list:
        protein_stats_map = {}
        progress_count += 1
        protein_mapped_peptides = []

        if protein.protein in proteins_to_peptide_map:
            protein_mapped_peptides = proteins_to_peptide_map[protein.protein]

        #Calculating other stats
        # print("Candidate Peptides Size", len(protein_mapped_peptides))
        coverage_array = protein.coverage_by_amino_acids(protein_mapped_peptides)
        coverage = float(sum(coverage_array))/float(len(coverage_array))

        #Checking out coverable amino acids
        uncoverable_array = protein.coverage_of_uncoverable_sequences()
        coverable_array = protein.coverage_of_coverable_sequences()
        uncoverable_size = sum(uncoverable_array)
        coverable_size = sum(coverable_array)

        covered_amino_acids_on_coverable = ming_numerical_utilities.dot_product(coverable_array, coverage_array)
        covered_amino_acids_on_uncoverable = ming_numerical_utilities.dot_product(uncoverable_array, coverage_array)
        coverage_of_coverable = 0.0
        coverage_of_uncoverable = 0.0
        if coverable_size > 0:
            coverage_of_coverable = float(covered_amino_acids_on_coverable)/float(coverable_size)
        if uncoverable_size > 0:
            coverage_of_uncoverable = float(covered_amino_acids_on_uncoverable)/float(uncoverable_size)
        coverable_percentage = float(coverable_size)/float(len(uncoverable_array))

        protein_stats_map["coverage_of_coverable"] = coverage_of_coverable
        protein_stats_map["aminoacid_coverage_of_coverable"] = covered_amino_acids_on_coverable
        protein_stats_map["coverage_of_uncoverable"] = coverage_of_uncoverable
        protein_stats_map["aminoacid_coverage_of_uncoverable"] = covered_amino_acids_on_uncoverable
        protein_stats_map["coverable_size"] = coverable_size
        protein_stats_map["uncoverable_size"] = uncoverable_size
        protein_stats_map["coverable_percentage"] = coverable_percentage


        protein_stats_map["aminoacid_coverage"] = coverage
        protein_stats_map["protein_length"] = len(coverage_array)
        # print(protein.protein + " " + str(coverage) + " " + str(progress_count) + " of " + str(len(protein_list)))

        number_of_peptides_covered = len(protein_mapped_peptides)
        protein_stats_map["peptides_on_protein"] = number_of_peptides_covered

        # print("number_of_peptides_covered", number_of_peptides_covered)
        try:
            uniprot_accession = protein.protein.split("|")[1]
        except:
            uniprot_accession = protein.protein
        # if uniprot_accession.find("-") != -1:
        protein_coverage_map[protein.protein] = protein_stats_map
        # continue
        protein_stats_map = uniprot_domain_statistics.calculate_protein_statistics(uniprot_accession, protein_to_annotations_map, protein, protein_mapped_peptides, coverage_array)

        protein_stats_map["coverage_of_coverable"] = coverage_of_coverable
        protein_stats_map["aminoacid_coverage_of_coverable"] = covered_amino_acids_on_coverable
        protein_stats_map["coverage_of_uncoverable"] = coverage_of_uncoverable
        protein_stats_map["aminoacid_coverage_of_uncoverable"] = covered_amino_acids_on_uncoverable
        protein_stats_map["coverable_size"] = coverable_size
        protein_stats_map["uncoverable_size"] = uncoverable_size
        protein_stats_map["coverable_percentage"] = coverable_percentage

        protein_stats_map["aminoacid_coverage"] = coverage
        protein_stats_map["peptides_on_protein"] = number_of_peptides_covered
        protein_stats_map["protein_length"] = len(coverage_array)
        protein_stats_map["covered_amino_acids"] = sum(coverage_array)
        protein_stats_map["pe_number"] = protein.pe_number
        protein_coverage_map[protein.protein] = protein_stats_map

        debug_coverage_dict[uniprot_accession] = coverage_array

    open(output_coverage_filename, "w").write(json.dumps(protein_coverage_map))

    #Outting debug
    output_debug_coverage_filename = os.path.join(output_protein_coverage_folder, str(partition_of_node) + ".json")
    open(output_debug_coverage_filename, "w").write(json.dumps(debug_coverage_dict))

    exit(0)


if __name__ == "__main__":
    main()
