#!/usr/bin/python


import sys
import getopt
import os
import json
import ming_fileio_library
import ming_protein_library
from collections import defaultdict

def usage():
    print("<input intermediate folder of protein coverage> <input intermediate folder of unique peptides> <output protein coverage file> <output unique peptides file> <library_coverage_summary_statistics>")


def hupo_nonoverlapping(segments, just_noncontained = True):
    if len(segments) == 0:
        return 0
    max_peps = [0 for _ in range(len(segments))]
    max_score_at_pep = [0 for _ in range(len(segments))]
    segments = [(int(s[0]),int(s[1])) for s in segments]
    ordered_segments = sorted(segments, key = lambda x: (int(x[1]),int(x[0])))

    for i in range(len(segments)):
        max_peps_j = 0
        max_score_at_pep_j = 0

        j = None
        for k in range(i-1,-1,-1):
            if (ordered_segments[k][0] <= ordered_segments[i][0] and
                ordered_segments[k][1] < ordered_segments[i][1] and
                ordered_segments[k][1] >= ordered_segments[i][0]):
                if j:
                    if (max_peps[k] > max_peps[j]):
                        j = k
                else:
                    j = k

        if j:
            max_peps_j = max_peps[j]

        j = None
        for k in range(i-1,-1,-1):
            if (ordered_segments[k][1] <= ordered_segments[i][0]):
                j = k
                break

        if j is not None:
            max_score_at_pep_j = max_score_at_pep[j]

        max_peps[i] = 1 + max(max_peps_j,max_score_at_pep_j)

        max_score_at_pep_j = 0

        j = None
        for k in range(i-1,-1,-1):
            if (ordered_segments[k][1] < ordered_segments[i][1] and
                ordered_segments[k][1] >= ordered_segments[i][0]):
                if j:
                    if (max_score_at_pep[k] > max_score_at_pep[j]):
                        j = k
                else:
                    j = k

        if j is not None:
            max_score_at_pep_j = max_score_at_pep[j]

        max_score_at_pep[i] = max(max_score_at_pep_j,max_peps[i])

    return max_score_at_pep[-1]


def summarize_statistics(proteome, peptides_to_protein_mapping, all_protein_stats, output_summary_filename):

    tryptic_peptides = proteome.get_tryptic_peptides()
    input_peptides = peptides_to_protein_mapping.keys()
    intersection_tryptic_peptides = list(set(tryptic_peptides).intersection(set(input_peptides)))

    total_covered_amino_acids = 0
    total_amino_acids_in_proteome = 0
    # for protein in all_protein_stats:
    #     total_covered_amino_acids += all_protein_stats[protein]["covered_amino_acids"]
    #     total_amino_acids_in_proteome += all_protein_stats[protein]["protein_length"]

    output_dict = defaultdict(list)

    output_dict["category"].append("Number unique sequnces")
    output_dict["value"].append(len(input_peptides))

    output_dict["category"].append("Number Tryptic sequnces")
    output_dict["value"].append(len(tryptic_peptides))

    output_dict["category"].append("Number Tryptic Sequences in Library")
    output_dict["value"].append(len(intersection_tryptic_peptides))

    output_dict["category"].append("Total Amino Acids Covered in Proteome")
    output_dict["value"].append(total_covered_amino_acids)

    output_dict["category"].append("Total Amino Acids in Proteome")
    output_dict["value"].append(proteome.get_proteome_amino_acid_count())

    ming_fileio_library.write_dictionary_table_data(output_dict, output_summary_filename)

def main():
    proteome_filename = sys.argv[1]
    input_intermediate_folder = sys.argv[2]
    input_intermediate_unique_peptides_folder = sys.argv[3]
    output_protein_coverage_filename = sys.argv[4]
    output_unique_peptides = sys.argv[5]
    input_nextprot_pe = sys.argv[6]
    library_coverage_summary_statistics = sys.argv[7]

    all_protein_stats = {}

    #Creating a command line for each partition
    all_intermediate_files = ming_fileio_library.list_files_in_dir(input_intermediate_folder)
    for parallel_output_filename in all_intermediate_files:
        print("Loading", parallel_output_filename)
        json_parition = json.loads(open(parallel_output_filename, "r").read())

        for key in json_parition:
            all_protein_stats[key] = json_parition[key]

    print "Total Proteins: " + str(len(all_protein_stats))


    #Loading the unique peptides per protein
    all_intermediate_peptide_mapping_files = ming_fileio_library.list_files_in_dir(input_intermediate_unique_peptides_folder)
    all_proteins_to_gene_unique_peptides_map = defaultdict(list)
    all_proteins_to_gene_unique_peptides_map_no_mismatch = defaultdict(list)
    all_proteins_to_protein_unique_peptides_map = defaultdict(list)
    all_proteins_to_protein_unique_peptides_map_no_mismatch = defaultdict(list)

    peptides_to_protein_mapping = defaultdict(list)

    for parallel_output_filename in all_intermediate_peptide_mapping_files:
        print("Loading", parallel_output_filename)
        json_parition = json.loads(open(parallel_output_filename, "r").read())

        for protein in json_parition["protein_to_gene_unique_peptide_no_mismatch"]:
            all_proteins_to_gene_unique_peptides_map_no_mismatch[protein] += json_parition["protein_to_gene_unique_peptide_no_mismatch"][protein]

        for protein in json_parition["protein_to_gene_unique_peptide"]:
            all_proteins_to_gene_unique_peptides_map[protein] += json_parition["protein_to_gene_unique_peptide"][protein]

        for protein in json_parition["protein_to_protein_unique_peptide_no_mismatch"]:
            all_proteins_to_protein_unique_peptides_map_no_mismatch[protein] += json_parition["protein_to_protein_unique_peptide_no_mismatch"][protein]

        for protein in json_parition["protein_to_protein_unique_peptide"]:
            all_proteins_to_protein_unique_peptides_map[protein] += json_parition["protein_to_protein_unique_peptide"][protein]

        for peptide in json_parition["peptide_to_protein_list"]:
            peptides_to_protein_mapping[peptide] = json_parition["peptide_to_protein_list"][peptide]

    for protein in all_protein_stats:
        if protein in all_proteins_to_gene_unique_peptides_map:
            all_protein_stats[protein]["gene_unique_peptide_count"] = min(len(set([p[0] for p in all_proteins_to_gene_unique_peptides_map[protein]])),hupo_nonoverlapping([p[1] for p in all_proteins_to_gene_unique_peptides_map[protein]]))
        else:
            all_protein_stats[protein]["gene_unique_peptide_count"] = 0
        if protein in all_proteins_to_gene_unique_peptides_map_no_mismatch:
            all_protein_stats[protein]["gene_unique_peptide_count_no_mismatch"] = len(set([p[0] for p in all_proteins_to_gene_unique_peptides_map_no_mismatch[protein]]))
        else:
            all_protein_stats[protein]["gene_unique_peptide_count_no_mismatch"] = 0

    for protein in all_protein_stats:
        if protein in all_proteins_to_protein_unique_peptides_map:
            all_protein_stats[protein]["unique_peptide_count"] = min(len(set([p[0] for p in all_proteins_to_protein_unique_peptides_map[protein]])),hupo_nonoverlapping([p[1] for p in all_proteins_to_protein_unique_peptides_map[protein]]))
        else:
            all_protein_stats[protein]["unique_peptide_count"] = 0
        if protein in all_proteins_to_protein_unique_peptides_map_no_mismatch:
            all_protein_stats[protein]["unique_peptide_count_no_mismatch"] = len(set([p[0] for p in all_proteins_to_protein_unique_peptides_map_no_mismatch[protein]]))
        else:
            all_protein_stats[protein]["unique_peptide_count_no_mismatch"] = 0

    #Format this appropriately
    output_dict = defaultdict(list)
    for protein in all_protein_stats:
        output_dict["protein"].append(protein)
        for key in all_protein_stats[protein]:
            output_dict[key].append(all_protein_stats[protein][key])
    ming_fileio_library.write_dictionary_table_data(output_dict, output_protein_coverage_filename)

    #Outputting all the unique peptides
    peptide_mapping_output_dict = defaultdict(list)
    for peptide in peptides_to_protein_mapping:
        peptide_mapping_output_dict["peptide"].append(peptide)
        peptide_mapping_output_dict["protein_count"].append(len(peptides_to_protein_mapping[peptide]))
        peptide_mapping_output_dict["proteins"].append(";".join(peptides_to_protein_mapping[peptide]))
    ming_fileio_library.write_dictionary_table_data(peptide_mapping_output_dict, output_unique_peptides)


    proteome = ming_protein_library.parse_fasta_proteome_file(proteome_filename,input_nextprot_pe)
    summarize_statistics(proteome, peptides_to_protein_mapping, all_protein_stats, library_coverage_summary_statistics)


if __name__ == "__main__":
    main()
