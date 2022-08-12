import sys
from collections import defaultdict
from pathlib import Path
import argparse
from python_ms_utilities import mapping, exon_mapping
from csv import DictWriter
from bisect import bisect_left

def arguments():
    parser = argparse.ArgumentParser(description='Map Peptides')
    parser.add_argument('-f','--proteome_fasta', type = Path, help='Input FASTA Protein Database')
    parser.add_argument('-c','--contaminants_fasta', type = Path, help='Input FASTA Protein Contaminants Database')
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

def main():
    args = arguments()

    proteome = mapping.merge_proteomes([mapping.read_uniprot(args.proteome_fasta),mapping.read_fasta(args.contaminants_fasta)])
    proteome_with_decoys = mapping.add_decoys(proteome)
    peptide_list = load_peptide_list(args.peptide_list)

    # peptide_to_exon_map = exon_mapping.map_peptides_to_exons(args.exon_fasta,peptide_list)
    peptide_to_exon_map = {}
    kmer_proteome_hashes = mapping.prepare_proteome_for_mapping(proteome_with_decoys, peptide_list)

    output_dict = defaultdict(list)

    fieldnames = [
        "peptide",
        "il_peptide",
        "mapped_proteins",
        "mapped_exons"
    ]

    with open(args.output_folder.joinpath(args.peptide_list.name), 'w') as f:
        w = DictWriter(f, delimiter = '\t', fieldnames = fieldnames)
        w.writeheader()

        for peptide in peptide_list:

            output_dict = {}

            mapped_proteins = mapping.map_peptide_to_proteome(peptide,proteome_with_decoys,kmer_proteome_hashes)
            mapped_exons = peptide_to_exon_map.get(peptide,[])

            output_dict["peptide"] = peptide
            output_dict["il_peptide"] = peptide.replace("I", "L")
            output_dict["mapped_proteins"] = mapping.protein_mappings_to_string(mapped_proteins)
            output_dict["mapped_exons"] = mapping.exon_mappings_to_string(mapped_exons)

            w.writerow(output_dict)

if __name__ == "__main__":
    main()
