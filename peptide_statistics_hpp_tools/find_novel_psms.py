import argparse
import sys
import csv
import glob
from quantlib import mztab
from collections import defaultdict


def arguments():
    parser = argparse.ArgumentParser(description='mzTab to list of peptides')
    parser.add_argument('-m','--mztab', type = str, help='mzTab')
    parser.add_argument('-s','--spec_on_server', type = str, help='Spectra')
    parser.add_argument('-n','--novel_coverage', type = str, help='Novel Coverage')
    parser.add_argument('-p','--novel_psms', type = str, help='Novel PSMs')
    if len(sys.argv) < 4:
        parser.print_help()
        sys.exit(1)
    return parser.parse_args()

def main():
    args = arguments()
    ids = {}
    if args.mztab != "" and args.mztab != None:
        for mztab_file in glob.glob(args.mztab + '/*'):
            ids = mztab.read(mztab_file, ids)
    peptide_to_psm = defaultdict(list)
    for filescan,psm in ids.items():
        peptide_to_psm[''.join([p.replace('I','L') for p in psm[0].sequence if p.isalpha()])].append((filescan,psm))

    with open(args.novel_coverage) as f, open(args.novel_psms, 'w') as fw:
        header = ['protein','filename','scan','sequence','charge','type']
        r = csv.DictReader(f, delimiter = '\t')
        w = csv.DictWriter(fw, delimiter = '\t', fieldnames = header)
        w.writeheader()
        for l in r:
            supporting_peptides = []
            novel_peptides = []
            if l['supporting_peptides'] != ' ':
                supporting_peptides = [p.split(' ')[0] for p in l['supporting_peptides'].split(',')]
            if l['novel_peptides'] != ' ':
                novel_peptides = [p.split(' ')[0] for p in l['novel_peptides'].split(',')]
            protein = l['protein']
            for peptide in supporting_peptides + novel_peptides:
                for (filescan,psm) in peptide_to_psm[peptide]:
                    w.writerow({
                        'protein':protein,
                        'filename':filescan[0],
                        'scan':filescan[1],
                        'sequence':psm[0].sequence,
                        'charge':psm[0].charge,
                        'type':'Supporting' if peptide in supporting_peptides else 'Novel'
                    })



if __name__ == '__main__':
    main()
