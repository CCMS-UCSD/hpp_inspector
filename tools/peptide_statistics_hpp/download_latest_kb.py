import argparse
import sys
import csv
from pathlib import Path
from collections import defaultdict
import read_mappings
from python_ms_utilities import mapping, resources
import pandas as pd 

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
    all_protein_mappings = []
    for protein_coverage_file in input_folder.glob('*'):
        _,protein_mappings_out,_ = read_mappings.read_protein_coverage(protein_coverage_file, proteome)
        all_protein_mappings.extend(protein_mappings_out)
    return pd.DataFrame(all_protein_mappings)

def main():
    args = arguments()

    if args.comparisons:

        try:

            proteome = mapping.read_uniprot(args.proteome_fasta)
            proteome = mapping.add_decoys(proteome)

            header = ['protein','aa_start','aa_end','demodified','synthetic_cosine']
            with open(args.kb_pep, 'w') as w:
                r = csv.DictWriter(w, delimiter = '\t', fieldnames = header)
                r.writeheader()
                for protein_mapping_obj in read_coverage_folder(args.comparisons, proteome).items():
                    r.writerow({
                        'protein':protein_mapping_obj['protein_accession'],
                        'aa_start':protein_mapping_obj['start_pos'],
                        'aa_end':protein_mapping_obj['end_pos'],
                        'demodified':protein_mapping_obj['peptide'],
                        'synthetic_cosine':protein_mapping_obj['synthetic_cosine'],
                    })
        except:
            header = ['protein','aa_start','aa_end','demodified','synthetic_cosine']
            with open(args.kb_pep, 'w') as w:
                r = csv.DictWriter(w, delimiter = '\t', fieldnames = header)
                r.writeheader()

                in_vivo = set()
                synthetic = set()

                for kb_input in args.comparisons.glob('*'):
                    with open(kb_input) as f:
                        kb_rs = csv.DictReader(f, delimiter = '\t')
                        for kb_row in kb_rs:
                            if kb_row['library'] == '2':
                                in_vivo.add((
                                    kb_row['protein'],
                                    kb_row['aa_start'],
                                    kb_row['aa_end'],
                                    kb_row['demodified']
                                ))
                            else:
                                synthetic.add((
                                    kb_row['protein'],
                                    kb_row['aa_start'],
                                    kb_row['aa_end'],
                                    kb_row['demodified']
                                ))
                for in_vivo_seq in in_vivo:
                    r.writerow({
                        'protein':in_vivo_seq[0],
                        'aa_start':in_vivo_seq[1],
                        'aa_end':in_vivo_seq[2],
                        'demodified':in_vivo_seq[3],
                        'synthetic_cosine': '0' if in_vivo_seq in synthetic else '-1'
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
