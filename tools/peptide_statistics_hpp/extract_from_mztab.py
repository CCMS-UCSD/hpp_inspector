import argparse
import sys
import csv
import glob
import requests
from collections import defaultdict
from quantlib import mztab
import urllib
import requests
import json
import time
from pathlib import Path
import pickle
import math

def arguments():
    parser = argparse.ArgumentParser(description='mzTab to list of peptides')
    parser.add_argument('-m','--mztab', type = Path, help='mzTab')
    parser.add_argument('-x','--novel_psms', type = Path, help='Novel PSMs')
    parser.add_argument('-p','--peptide_folder', type = Path, help='Output Peptide Folder')
    parser.add_argument('-l','--mapping_parallelism', type = int, help='Parallelism for Mapping')
    parser.add_argument('-t','--peak_tolerance', type = float, help='Peak Tolerance for Matching')

    if len(sys.argv) < 3:
        parser.print_help()
        sys.exit(1)
    return parser.parse_args()

def main():
    args = arguments()

    ids = {}

    if args.mztab != None:
        for mztab_file in args.mztab.glob('*'):
            try:
                ids = mztab.read(mztab_file, ids)
            except:
                ids = mztab.read_lib(mztab_file, ids)

    peptides = set()

    with open(args.novel_psms, 'w') as fw_psm:
        header = ['filename','scan','sequence','charge','score','pass','type','parent_mass','frag_tol']
        w_psm = csv.DictWriter(fw_psm, delimiter = '\t', fieldnames = header)
        w_psm.writeheader()
        for (filescan,psm) in ids.items():
            peptides.add(''.join([aa for aa in psm[0].sequence if aa.isalpha()]))
            psm_row = {
                'filename':filescan[0],
                'scan':filescan[1],
                'sequence':psm[0].sequence,
                'charge':psm[0].charge,
                'score':psm[0].score,
                'pass':'N/A',
                'parent_mass':psm[0].parent_mass,
                'frag_tol':psm[0].fragment_tolerance if psm[0].fragment_tolerance else float(args.peak_tolerance)
            }
            w_psm.writerow(psm_row)

    peptides = list(peptides)
    if len(peptides) < args.mapping_parallelism:
        n = args.mapping_parallelism
    else:
        n = math.ceil(len(peptides)/args.mapping_parallelism)
    split_peptides = [peptides[i:i+n] for i in range(0, len(peptides), n)]

    for i, split_peptide_list in enumerate(split_peptides):
        with open(args.peptide_folder.joinpath(str(i)), 'w') as w:
            for peptide in split_peptide_list:
                w.write("{}\n".format(peptide))

if __name__ == '__main__':
    main()
