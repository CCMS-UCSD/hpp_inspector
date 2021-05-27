import argparse
import sys
import csv
from pathlib import Path
from collections import defaultdict
import read_mappings

def arguments():
    parser = argparse.ArgumentParser(description='mzTab to list of peptides')
    parser.add_argument('-p','--params', type = str, help='Input Parameters')
    parser.add_argument('-c','--comparisons', type = Path, help='Comparison Jobs')
    parser.add_argument('-b','--backup_kb_pep', type = Path, help='Backup KB Peptides')
    parser.add_argument('-k','--kb_pep', type = Path, help='Output KB Peptides')
    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit(1)
    return parser.parse_args()

def read_coverage_folder(input_folder):
    nextprot_pe = {}
    added_proteins = defaultdict(lambda: defaultdict(list))
    all_proteins = defaultdict(lambda: defaultdict(list))
    pep_mapping_info = {}
    peptide_to_exon_map = defaultdict(list)
    for protein_coverage_file in input_folder.glob('*'):
        all_proteins,added_proteins,pep_mapping_info,peptide_to_exon_map = read_mappings.read_protein_coverage(protein_coverage_file, all_proteins,added_proteins,pep_mapping_info,peptide_to_exon_map,nextprot_pe)
    return added_proteins

def main():
    args = arguments()

    if args.comparisons:
        header = ['library','protein','aa_start','aa_end','demodified']
        with open(args.kb_pep, 'w') as w:
            r = csv.DictWriter(w, delimiter = '\t', fieldnames = header)
            r.writeheader()
            for protein, peptide_mappings in read_coverage_folder(args.comparisons).items():
                print(protein)
                print(peptide_mappings)
                for peptide, mappings in peptide_mappings.items():
                    for (start, end) in mappings:
                        r.writerow({
                            'library':2,
                            'protein':protein,
                            'aa_start':start,
                            'aa_end':end,
                            'demodified':peptide
                        })
    else:
        args.kb_pep.write_text(args.backup_kb_pep.read_text())

    # try:
    #
    #     PAGE_SIZE = 500000;
    #     OFFSET = 0;
    #
    #     peptides_to_add = requests.get('http://ccms-internal.ucsd.edu/ProteoSAFe/PeptideStatisticsServlet?pageSize={}&offset={}'.format(PAGE_SIZE,OFFSET)).json()
    #     peptides = peptides_to_add
    #     while (len(peptides_to_add) == PAGE_SIZE):
    #         OFFSET += PAGE_SIZE;
    #         peptides_to_add = requests.get('http://ccms-internal.ucsd.edu/ProteoSAFe/PeptideStatisticsServlet?pageSize={}&offset={}'.format(PAGE_SIZE,OFFSET)).json()
    #         peptides.extend(peptides_to_add)
    #
    #     header = ['library','protein','aa_start','aa_end','demodified','charge', 'sequence']
    #     with open(args.kb_pep, 'w') as w:
    #         r = csv.DictWriter(w, delimiter = '\t', fieldnames = header)
    #         r.writeheader()
    #         for peptide in peptides:
    #             r.writerow(peptide)
    #
    # except:



if __name__ == '__main__':
    main()
