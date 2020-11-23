#!/usr/bin/python


import sys
import getopt
import os
import json
import ming_protein_library
import ming_fileio_library
import re
import uniprot_domain_statistics
from collections import defaultdict

if (sys.version_info > (3, 0)):
    # Python 3 code in this block
    import _pickle as pickle
else:
    # Python 2 code in this block
    import cPickle as pickle

def usage():
    print("<input fasta> <input peptides> <input parameters file> <path to annotations> <output results folder>")

def strip_sequence(sequence):
    if sequence[1] == "." and sequence[-2] == ".":
        sequence = sequence[2:-2]
    # whitelist =  ["57.021464","15.994915","14.015650","17.026549","43.005814","42.010565","0.984016"]
    # whitelist += ["57.021",   "15.995",   "14.016",   "17.027",   "43.006",   "42.011",   "0.984"   ,"229.163"]
    #
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
        return_peptides.append(peptide)


    return_peptides = list(set(return_peptides))
    return_peptides.sort()

    return return_peptides

def main():
    input_fasta_filename = sys.argv[1]
    input_peptide_list_filename = sys.argv[2]
    input_parameters_filename = sys.argv[3]
    input_uniprot_gff_pkl_filename = sys.argv[4]
    input_nextprot_pe = sys.argv[5]
    output_results_folder = sys.argv[6]

    partition_total = 1
    partition_of_node = 0

    print("LOADING " + input_parameters_filename)

    params_map = json.loads(open(input_parameters_filename).read())
    partition_total = params_map["total_paritions"]
    partition_of_node = params_map["node_partition"]

    proteome = ming_protein_library.parse_fasta_proteome_file(input_fasta_filename, input_nextprot_pe)
    #IL substitution for proteome

    protein_to_gene = {}

    for protein in proteome.protein_list:
        protein.sequence = protein.sequence.replace("I", "L")
        protein_to_gene[protein.protein] = protein.gene_name

    peptide_list = load_peptide_list(input_peptide_list_filename)
    protein_list = proteome.protein_list

    #Lets alternate what we get
    peptide_list = list(set([strip_sequence(x) for x in peptide_list[partition_of_node::partition_total] if strip_sequence(x) != None]))

    #IL peptides
    il_peptide_list = list(set([x.replace("I", "L") for x in peptide_list]))

    protein_to_unique_peptide_map = defaultdict(list)
    peptide_to_protein_map_w_mismatch = proteome.get_peptides_mapped_to_proteins_mismatched(il_peptide_list, mismatch_max=1)
    peptide_to_protein_map = proteome.get_peptides_mapped_to_proteins_efficient_w_coordinates(il_peptide_list)

    print("Loading PKL of Uniprot annotations")
    protein_to_annotations_map = pickle.load(open(input_uniprot_gff_pkl_filename, "rb"))

    output_dict = defaultdict(list)

    additional_columns = ["protein_name", "gene", "activesites_covered", "mod_covered", "metalbinding_covered", "glycosylation_covered", "bindingsites_covered", "disulfide_covered"]

    count = 0

    print(len(peptide_list))
    print(len(set(peptide_list)))

    for peptide in peptide_list:

        if strip_sequence(peptide):
            count += 1

            if count % 100 == 0:
                print(count, "of", len(peptide_list))


            proteins_covered_w_coords = peptide_to_protein_map[strip_sequence(peptide).replace("I", "L")]
            proteins_covered_w_coords_mismatch = peptide_to_protein_map_w_mismatch[strip_sequence(peptide).replace("I", "L")]
            proteins_covered = list(set([p.split(" (")[0] for p in proteins_covered_w_coords]))
            proteins_covered_all = list(set([p.split(" (")[0] for p in proteins_covered_w_coords_mismatch]))
            num_proteins_covered = len(set(peptide_to_protein_map_w_mismatch[strip_sequence(peptide).replace("I", "L")]))
            proteins_covered_no_iso_no_tr = set([p.split("|")[1].split("-")[0] for p in proteins_covered_all if "sp" in p])
            proteins_covered_no_tr = set([p for p in proteins_covered_all if "sp" in p])
            output_dict["peptide"].append(peptide)
            output_dict["il_peptide"].append(strip_sequence(peptide).replace("I", "L"))
            output_dict["proteins_w_coords"].append(";".join(proteins_covered_w_coords))
            output_dict["proteins"].append(";".join(proteins_covered))
            output_dict["proteins_all"].append(";".join(proteins_covered_all))
            output_dict["proteins_no_iso_no_tr"].append(";".join(proteins_covered_no_iso_no_tr))
            output_dict["proteins_no_tr"].append(";".join(proteins_covered_no_tr))
            output_dict["proteinsg_covered"].append(set([protein_to_gene[p] for p in proteins_covered]))
            output_dict["proteinsg_covered_mismatch"].append(set([protein_to_gene[p] for p in proteins_covered_all]))
            output_dict["num_proteins_all"].append(len(proteins_covered_all))
            output_dict["num_proteins_no_iso_no_tr"].append(len(proteins_covered_no_iso_no_tr))
            output_dict["num_proteins_no_tr"].append(len(proteins_covered_no_tr))
            output_dict["num_proteins"].append(len(proteins_covered))
            output_dict["num_proteins_w_mismatch"].append(num_proteins_covered)
            output_dict["num_proteinsg_covered"].append(len(set([protein_to_gene[p] for p in proteins_covered])))

            shared_gene = True
            min_genes_covered = None
            for protein_pos in peptide_to_protein_map[strip_sequence(peptide).replace("I", "L")]:
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
            output_dict["gene_unique_no_mismatch"].append(shared_gene)
            for p in peptide_to_protein_map_w_mismatch[strip_sequence(peptide).replace("I", "L")]:
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
            output_dict["gene_unique"].append(shared_gene)
            pe_covered = [proteome.protein_map[protein.split(" ")[0]].pe_number for protein in proteins_covered]
            output_dict["pe1"].append(sum([1 if pe == "1" else 0 for pe in pe_covered]))
            output_dict["pe2"].append(sum([1 if pe == "2" else 0 for pe in pe_covered]))
            output_dict["pe3"].append(sum([1 if pe == "3" else 0 for pe in pe_covered]))
            output_dict["pe4"].append(sum([1 if pe == "4" else 0 for pe in pe_covered]))
            output_dict["pe5"].append(sum([1 if pe == "5" else 0 for pe in pe_covered]))

            if len(set(peptide_to_protein_map[strip_sequence(peptide).replace("I", "L")])) != 1:
                for key in additional_columns:
                    if key == 'protein_name':
                        output_dict[key].append("Multiple proteins")
                    elif key == 'gene':
                        output_dict[key].append("Multiple genes")
                    else:
                        output_dict[key].append(0)
                continue

            protein_object = proteome.protein_map[proteins_covered[0].split(" ")[0]]
            place_holder_peptide_list = [strip_sequence(peptide).replace("I", "L")]
            coverage_array = protein_object.coverage_by_amino_acids(place_holder_peptide_list)
            try:
                uniprot_accession = protein_object.protein.split("|")[1]
            except:
                uniprot_accession = protein_object.protein
            # if uniprot_accession.find("-") != -1:
            #     for key in additional_columns:
            #         output_dict[key].append(0)
            #     continue

            protein_stats_map = uniprot_domain_statistics.calculate_protein_statistics(uniprot_accession, protein_to_annotations_map, protein_object, place_holder_peptide_list, coverage_array)
            for key in additional_columns:
                output_dict[key].append(protein_stats_map[key])



    output_coverage_filename = os.path.join(output_results_folder, str(partition_of_node) + ".tsv")
    ming_fileio_library.write_dictionary_table_data(output_dict, output_coverage_filename)


if __name__ == "__main__":
    main()
