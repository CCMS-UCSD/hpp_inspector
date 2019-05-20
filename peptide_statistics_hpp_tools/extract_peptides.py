import argparse
import sys
import csv
import glob
import mztab

def arguments():
    parser = argparse.ArgumentParser(description='mzTab to list of peptides')
    parser.add_argument('-i','--mztab', type = str, help='Input Identifications')
    parser.add_argument('-p','--peptidelist', type = str, help='Output Peptides')
    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit(1)
    return parser.parse_args()

def main():
    args = arguments()
    ids = {}
    if args.mztab != "" and args.mztab != None:
        for mztab_file in glob.glob(args.mztab + '/*'):
            ids = mztab.read(mztab_file, ids)
    peptides = set([
        ''.join([aa for aa in id[0].sequence if aa.isalpha()])
        for id in
        ids.values()
    ])
    with open(args.peptidelist,'w') as w:
        w.write('Peptide\n')
        for peptide in peptides:
            w.write('{}\n'.format(peptide))

if __name__ == '__main__':
    main()
