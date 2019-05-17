#!/usr/bin/python


import sys
import getopt
import os
import ming_protein_library
import re

def usage():
    print "<input fasta> <input peptide list> <output statistics>"

def main():
    input_fasta_filename = sys.argv[1]
    input_peptide_list = sys.argv[2]
    input_nextprot_pe = sys.argv[3]
    output_file = open(sys.argv[4], "w")

    min_length = 0

    proteome = ming_protein_library.parse_fasta_proteome_file(input_fasta_filename,input_nextprot_pe)
    peptide_list = []

    for line in open(input_peptide_list, "r"):
        peptide = line.rstrip()
        #peptide_list += re.findall(r".(?:(?<![KR](?!P)).)*", peptide)
        if len(peptide) > min_length:
            peptide_list.append(peptide)
    peptide_list = list(set(peptide_list))

    #Getting all trypic peptides for proteome
    tryptic_peptides = proteome.get_tryptic_peptides(min_length)

    output_string = ""

    output_string += "Number of Unique Input Sequences: " + str(len(peptide_list)) + "\n"

    output_string += "Number of Tryptic Peptides in Proteome: " + str(len(tryptic_peptides)) + "\n"

    proteome_covered, library_covered = proteome.calculate_tryptic_peptides_covered(peptide_list)
    output_string += "Coverage of all trypic coverage in proteome by input peptides: " + str(proteome_covered) + "\n"


    output_string += "Coverage of all library peptides by proteome: " + str(library_covered) + "\n"

    #intersection_set, proteome_exclusive, cover_exclusive = proteome.set_coverage_with_proteome(peptide_list)
    #open("cover_exclusive.out", "w").write(str(cover_exclusive))

    output_string += "Number of Non Shared Peptides in Proteome: " + str(len(proteome.get_unique_peptides(min_length))) + "\n"

    output_string += "Coverage of Unique Tryptic Peptides in Proteome: "  + str(proteome.calculate_unique_peptide_coverage(peptide_list)) + "\n"

    output_file.write(output_string)

if __name__ == "__main__":
    main()
