import argparse
import sys
import csv
from pathlib import Path
from collections import defaultdict
import read_mappings
from python_ms_utilities import mapping, resources

def arguments():
    parser = argparse.ArgumentParser(description='mzTab to list of peptides')
    parser.add_argument('-p','--params', type = str, help='Input Parameters')
    parser.add_argument('-c','--comparisons', type = Path, help='Comparison Jobs')
    parser.add_argument('-f','--proteome_fasta', type = Path, help='FASTA File')
    parser.add_argument('-b','--backup_kb_pep', type = Path, help='Backup KB Peptides')
    parser.add_argument('-k','--kb_pep', type = Path, help='Output KB Peptides')
    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit(1)
    return parser.parse_args()

def read_coverage_folder(input_folder,proteome):
    added_proteins = defaultdict(lambda: defaultdict(list))
    all_proteins = defaultdict(lambda: defaultdict(list))
    pep_mapping_info = {}
    peptide_to_exon_map = defaultdict(list)
    for protein_coverage_file in input_folder.glob('*'):
        all_proteins,added_proteins,pep_mapping_info,peptide_to_exon_map = read_mappings.read_protein_coverage(protein_coverage_file, all_proteins,added_proteins,pep_mapping_info,peptide_to_exon_map,proteome)
    return added_proteins

def main():
    args = arguments()

    if args.comparisons:

        try:

            proteome = mapping.read_uniprot(args.proteome_fasta)
            proteome = mapping.add_decoys(proteome)

            header = ['library','protein','aa_start','aa_end','demodified']
            with open(args.kb_pep, 'w') as w:
                r = csv.DictWriter(w, delimiter = '\t', fieldnames = header)
                r.writeheader()
                for protein, peptide_mappings in read_coverage_folder(args.comparisons, proteome).items():
                    print(protein)
                    print(peptide_mappings)
                    for peptide, mappings in peptide_mappings.items():
                        for (start, end) in mappings:
                            r.writerow({
                                'protein':protein,
                                'aa_start':start,
                                'aa_end':end,
                                'demodified':peptide
                            })
        except:
            header = ['library','protein','aa_start','aa_end','demodified']
            with open(args.kb_pep, 'w') as w:
                r = csv.DictWriter(w, delimiter = '\t', fieldnames = header)
                r.writeheader()
                for kb_input in args.comparisons.glob('*'):
                    with open(kb_input) as f:
                        kb_rs = csv.DictReader(f)
                        for kb_row in kb_rs:
                            if kb_row['library'] == '2':
                                r.writerow({
                                    'protein':kb_row['protein'],
                                    'aa_start':kb_row['aa_start'],
                                    'aa_end':kb_row['aa_end'],
                                    'demodified':kb_row['demodified']
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
