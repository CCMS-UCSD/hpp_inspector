import argparse
import sys
import csv
import glob
import requests
from collections import defaultdict
from quantlib import mztab
import urllib
import requests
import spectrum_alignment as sa
import json
import time
from pyteomics import mzml
import pathlib
import pickle


CUTOFF = {"6": 2.12167006027, "7": 2.32154253572, "8": 2.53031824698, "9": 2.7717776307, "10": 3.0419577266, "11": 3.20404758311, "12": 3.4028353721, "13": 3.14854063781, "14": 3.47053358426, "15": 3.35416414915, "16": 3.34418455532, "17": 3.24904542906, "18": 3.18426919317, "19": 3.01362943461, "20": 3.12316407632, "21": 3.08160158158, "22": 2.59466460406, "23": 2.96230440256, "24": 3.11610830789, "25": 2.93990420679, "26": 2.6947192629, "27": 2.43042157531, "28": 2.5287628139300002, "29": 2.26034401643, "30": 2.60979068254, "31": 2.70004841417, "32": 2.69919517925, "33": 2.18110553715, "34": 1.90115111418, "35": 1.5402648112000001, "36": 1.74919562203, "37": 1.88066887473, "38": 1.58471929702, "39": 1.73377627878, "40": 3.10312899149}

def arguments():
    parser = argparse.ArgumentParser(description='mzTab to list of peptides')
    parser.add_argument('-m','--mztab', type = str, help='mzTab')
    parser.add_argument('-s','--spec_on_server', type = str, help='Spectra')
    parser.add_argument('-x','--novel_psms', type = str, help='Novel PSMs')
    parser.add_argument('-e','--novel_peptides', type = str, help='Novel Peptides')
    parser.add_argument('-y','--synthetics', type = str, help='Matched Synthetics')
    parser.add_argument('-t','--peak_tolerance', type = str, help='Peak Tolerance for Matching')

    if len(sys.argv) < 4:
        parser.print_help()
        sys.exit(1)
    return parser.parse_args()

def main():
    args = arguments()
    nextprot_pe = defaultdict(lambda: 0)
    protein_length = defaultdict(lambda: 0)

    tol = float(args.peak_tolerance)

    ids = {}

    if args.mztab != "" and args.mztab != None:
        for mztab_file in glob.glob(args.mztab + '/*'):
            try:
                ids = mztab.read(mztab_file, ids)
            except:
                ids = mztab.read_lib(mztab_file, ids)

    peptide_to_psm = defaultdict(list)
    for filescan,psm in ids.items():
        #Don't include improper bioplex files - seems to be fixed
        # if not ('j10022_RNF14' in filescan[0] or 'j12271_TMPRSS3' in filescan[0] or 'j13961_ADAM32' in filescan[0] or 'j9837_ADCYAP1' in filescan[0] or 'q8920_ACCN4' in filescan[0]):
        peptide_to_psm[''.join([p.replace('I','L') for p in psm[0].sequence if p.isalpha()])].append((filescan,psm))

    supporting_peptides = set()

    psms_to_consider = defaultdict(set)

    with open(args.novel_psms, 'w') as fw_psm:
        header = ['filename','scan','sequence','charge','score','pass','type','parent_mass','frag_tol']
        w_psm = csv.DictWriter(fw_psm, delimiter = '\t', fieldnames = header)
        w_psm.writeheader()
        for peptide, psms in peptide_to_psm.items():
            for (filescan,psm) in psms:
                psm_row = {
                    'filename':filescan[0],
                    'scan':filescan[1],
                    'sequence':psm[0].sequence,
                    'charge':psm[0].charge,
                    'score':psm[0].score,
                    'pass':'Above' if psm[0].score > CUTOFF.get(str(len(peptide)),0) else 'Below',
                    'parent_mass':psm[0].parent_mass,
                    'frag_tol':psm[0].fragment_tolerance
                }
                w_psm.writerow(psm_row)

if __name__ == '__main__':
    main()
