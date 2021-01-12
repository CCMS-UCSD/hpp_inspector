#!/usr/bin/python

import sys
import getopt
import os
import json
import ming_protein_library
import ming_fileio_library
import re
from collections import defaultdict
from pathlib import Path
import argparse
import exon_mapping
from csv import DictReader, DictWriter
import math
from bisect import bisect_left

def arguments():
    parser = argparse.ArgumentParser(description='Map Peptides')
    parser.add_argument('-f','--proteome_fasta', type = Path, help='Input FASTA Protein Database')
    parser.add_argument('-e','--exon_fasta', type = Path, help='Input FASTA Exon Mapping Database')
    parser.add_argument('-p','--peptide_list', type = Path, help='Peptide List')
    parser.add_argument('-o','--output_folder', type = Path, help='Output Folder')

    if len(sys.argv) < 3:
        parser.print_help()
        sys.exit(1)
    return parser.parse_args()

def load_peptide_list(input_peptide_filename):
    peptides = []
    with open(input_peptide_filename) as f:
        for l in f:
            peptides.append(l.rstrip())
    return peptides

def exon_str(exon):
    chr, mapped, complete = exon
    mapped_str = '-'.join(["{}/{}".format(m[0],m[1]) for m in mapped])
    complete_str = '-'.join(["{}/{}".format(c[0],c[1]) for c in complete])
    return "{0}:{1}:{2}".format(chr, complete_str, mapped_str)

def main():
    args = arguments()

    #map to all protiens
    proteome = ming_protein_library.parse_fasta_proteome_file(args.proteome_fasta)
    proteome.create_decoy_proteins()

    #IL substitution for proteome
    protein_to_gene = {}

    for protein in proteome.protein_list:
        protein_to_gene[protein.protein] = protein.gene_name

    peptide_list = load_peptide_list(args.peptide_list)

    peptide_to_protein_map = proteome.get_peptides_mapped_to_proteins_efficient_w_coordinates_w_mismatch(peptide_list)
    peptide_to_exon_map = exon_mapping.map_peptides_to_exons(args.exon_fasta,peptide_list)

    output_dict = defaultdict(list)

    fieldnames = [
        "peptide",
        "il_peptide",
        "proteins",
        "proteins_w_coords",
        "proteins_all",
        "proteins_no_iso_no_tr",
        "proteins_no_tr",
        "proteinsg_covered",
        "proteinsg_covered_mismatch",
        "num_proteins_all",
        "num_proteins_no_iso_no_tr",
        "num_proteins_no_tr",
        "num_proteins",
        "num_proteins_w_mismatch",
        "num_proteinsg_covered",
        "gene_unique",
        "gene_unique_no_mismatch",
        "mapped_exons",
        "exon_junctions_covered",
        "exons_covered_no_junction",
        "total_unique_exons_covered"
    ]

    count = 0

    mapping_str = lambda m: "{} ({}-{}){}".format(m.protein,m.start_pos,m.end_pos, 'A' if m.il_ambiguous else '')
    # exon_str = lambda m: "{}:{}/{} ({}-{})".format(m[0], m[1][0], m[1][1], m[2][0], m[2][1])

    with open(args.output_folder.joinpath(args.peptide_list.name), 'w') as f:
        w = DictWriter(f, delimiter = '\t', fieldnames = fieldnames)
        w.writeheader()

        for peptide in peptide_list:
            oudput_dict = {}

            # proteome mapping

            complete_proteins_covered = peptide_to_protein_map[peptide]

            proteins_covered_exact = list(set([p.protein for p in complete_proteins_covered if p.mismatches == 0]))
            proteins_covered_mismatch = list(set([p.protein for p in complete_proteins_covered if p.mismatches == 1]))
            proteins_covered_no_iso_no_tr = set([p.split("|")[1].split("-")[0] for p in proteins_covered_mismatch if "sp" in p])
            proteins_covered_no_tr = set([p for p in proteins_covered_mismatch if "sp" in p])

            output_dict["peptide"] = peptide
            output_dict["il_peptide"] = peptide.replace("I", "L")
            output_dict["proteins_w_coords"] = ";".join([mapping_str(p) for p in complete_proteins_covered if p.mismatches == 0])
            output_dict["proteins"] = ";".join(proteins_covered_exact)
            output_dict["proteins_all"] = ";".join(proteins_covered_mismatch)
            output_dict["proteins_no_iso_no_tr"] = ";".join(proteins_covered_no_iso_no_tr)
            output_dict["proteins_no_tr"] = ";".join(proteins_covered_no_tr)
            output_dict["proteinsg_covered"] = set([protein_to_gene[p] for p in proteins_covered_exact])
            output_dict["proteinsg_covered_mismatch"] = set([protein_to_gene[p] for p in proteins_covered_mismatch])
            output_dict["num_proteins_all"] = len(proteins_covered_mismatch)
            output_dict["num_proteins_no_iso_no_tr"] = len(proteins_covered_no_iso_no_tr)
            output_dict["num_proteins_no_tr"] = len(proteins_covered_no_tr)
            output_dict["num_proteins"] = len(proteins_covered_exact)
            output_dict["num_proteins_w_mismatch"] = len(proteins_covered_mismatch)
            output_dict["num_proteinsg_covered"] = len(set([protein_to_gene[p] for p in proteins_covered_exact]))

            shared_gene = True
            min_genes_covered = None
            for p in proteins_covered_exact:
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
            output_dict["gene_unique_no_mismatch"] = shared_gene

            for p in proteins_covered_mismatch:
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
            output_dict["gene_unique"] = shared_gene

            mapped_exons = peptide_to_exon_map[peptide]

            output_dict["mapped_exons"] = ";".join([exon_str(exon) for exon in mapped_exons])
            output_dict["exon_junctions_covered"] = len([exon for exon in mapped_exons if len(exon[1]) > 1])
            output_dict["exons_covered_no_junction"] = len([exon for exon in mapped_exons if len(exon[1]) == 1])
            output_dict["total_unique_exons_covered"] = len(set([e for exon in mapped_exons for e in exon[2]]))
    # for peptide, mapped_exons in peptide_to_exon_map.items():
    #     print(peptide)
    #     for exon, mapping_locations in mapped_exons.items():
    #         print("\t{0}:{1}/{2}".format(*exon))
    #         for location in mapping_locations:
    #             print("\t\t{0} {1} ({2[0]}/{2[1]})".format(*location))
            # # exon mapping
            # potential_locations = set()
            # exons_mapped = peptide_to_exon_map[peptide]
            # print(exons_mapped)
            # for exon in exons_mapped:
            #     if exon.mismatches == 0:
            #         psplit = exon.protein.split("@")
            #         chromosome = psplit[0]
            #         transcript = psplit[1]
            #         gene = psplit[2]
            #         print(exon_map_select[gene][transcript])
            #         for location in exon_map_select[gene][transcript]:
            #             potential_locations.add((gene,location))
            #
            # exons_covered = set()
            # exons_fully_covered = 0
            # exons_partial_lengths = set()
            # exons_junctions_covered = 0
            #
            # exon_locations_covered = peptides_to_exons[peptide]
            # for (chromosome, (start,end)) in exon_locations_covered:
            #
            #     for (gene,(start_ref,end_ref)) in potential_locations:
            #         if start >= start_ref and end <= end_ref:
            #             if (end - start) != len(peptide) * 3:
            #                 exons_partial_lengths.add((start,end))
            #             exons_covered.add((chromosome, (start_ref,end_ref),(start,end)))
            # if len(exons_partial_lengths) > 0:
            #     exons_junctions_covered = math.ceil(len(exons_partial_lengths) / sum([e[1]-e[0] for e in exons_partial_lengths])/(len(peptide) * 3))
            #
            # output_dict["exons_mapped"] = ";".join([e.protein.replace("GFFtoFASTA@","") for e in exons_mapped if e.mismatches == 0])
            # output_dict["num_exons_mapped"] = len(exons_mapped)
            # output_dict["exons"] = ";".join([exon_str(exon) for exon in exons_covered])
            # output_dict["exons_covered"] = len(exons_covered)
            # output_dict["exon_junctions_covered"] = exons_junctions_covered

            w.writerow(output_dict)

if __name__ == "__main__":
    main()
