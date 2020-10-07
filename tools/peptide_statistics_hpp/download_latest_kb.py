import argparse
import sys
import csv
import glob
from quantlib import mztab
import requests

def arguments():
    parser = argparse.ArgumentParser(description='mzTab to list of peptides')
    parser.add_argument('-p','--params', type = str, help='Input Parameters')
    parser.add_argument('-k','--kb_pep', type = str, help='Output KB Peptides')
    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit(1)
    return parser.parse_args()

def main():
    args = arguments()

    PAGE_SIZE = 500000;
    OFFSET = 0;

    peptides_to_add = requests.get('http://ccms-internal.ucsd.edu/ProteoSAFe/PeptideStatisticsServlet?pageSize={}&offset={}'.format(PAGE_SIZE,OFFSET)).json()
    peptides = peptides_to_add
    while (len(peptides_to_add) == PAGE_SIZE):
        OFFSET += PAGE_SIZE;
        peptides_to_add = requests.get('http://ccms-internal.ucsd.edu/ProteoSAFe/PeptideStatisticsServlet?pageSize={}&offset={}'.format(PAGE_SIZE,OFFSET)).json()
        peptides.extend(peptides_to_add)

    header = ['library','protein','aa_start','aa_end','demodified','charge', 'sequence']
    with open(args.kb_pep, 'w') as w:
        r = csv.DictWriter(w, delimiter = '\t', fieldnames = header)
        r.writeheader()
        for peptide in peptides:
            r.writerow(peptide)


if __name__ == '__main__':
    main()
