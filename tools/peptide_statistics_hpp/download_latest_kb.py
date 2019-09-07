import argparse
import sys
import csv
import glob
import mztab
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
    peptides = requests.get('http://ccms-internal.ucsd.edu/ProteoSAFe/PeptideStatisticsServlet').json()
    header = ['library','protein','aa_start','aa_end','demodified','charge', 'sequence']
    with open(args.kb_pep, 'w') as w:
        r = csv.DictWriter(w, delimiter = '\t', fieldnames = header)
        r.writeheader()
        for peptide in peptides:
            r.writerow(peptide)


if __name__ == '__main__':
    main()
