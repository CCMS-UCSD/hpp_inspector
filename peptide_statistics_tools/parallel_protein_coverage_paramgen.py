#!/usr/bin/python


import sys
import getopt
import os
import json

def usage():
    print("<input fasta> <input peptides> <output parameters files> <parallelism>")


def load_peptide_list(input_peptide_filename):
    all_peptides = []
    for line in open(input_peptide_filename):
        peptide = line.rstrip()
        all_peptides.append(peptide)
    return all_peptides

def map_sequence_to_proteins(input_list):
    return input_list[0].get_proteins_with_sequence(input_list[1])

def main():
    print sys.argv

    input_fasta_filename = sys.argv[1]
    input_peptide_list_filename = sys.argv[2]
    output_parameters_folder = sys.argv[3]

    partition_count = int(sys.argv[4])

    #Creating a command line for each partition
    for i in range(partition_count):
        output_param_filename = os.path.join(output_parameters_folder, str(i) + ".json")
        output_map = {}
        output_map["total_paritions"] = partition_count
        output_map["node_partition"] = i
        open(output_param_filename, "w").write(json.dumps(output_map))

    exit(0)


if __name__ == "__main__":
    main()
