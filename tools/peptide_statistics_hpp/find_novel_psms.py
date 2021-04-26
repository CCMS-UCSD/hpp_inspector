import argparse
import sys
import csv
import glob
from quantlib import mztab
from collections import defaultdict
import urllib
import requests
import spectrum_alignment as sa
import json
import time
from pyteomics import mzml
import pathlib
import pickle

CUTOFF = {"6": 2.12167006027, "7": 2.32154253572, "8": 2.53031824698, "9": 2.7717776307, "10": 3.0419577266, "11": 3.20404758311, "12": 3.4028353721, "13": 3.14854063781, "14": 3.47053358426, "15": 3.35416414915, "16": 3.34418455532, "17": 3.24904542906, "18": 3.18426919317, "19": 3.01362943461, "20": 3.12316407632, "21": 3.08160158158, "22": 2.59466460406, "23": 2.96230440256, "24": 3.11610830789, "25": 2.93990420679, "26": 2.6947192629, "27": 2.43042157531, "28": 2.5287628139300002, "29": 2.26034401643, "30": 2.60979068254, "31": 2.70004841417, "32": 2.69919517925, "33": 2.18110553715, "34": 1.90115111418, "35": 1.5402648112000001, "36": 1.74919562203, "37": 1.88066887473, "38": 1.58471929702, "39": 1.73377627878, "40": 3.10312899149}

# def get_peaks(filename, scan):
#     url = create_peaks_url(filename,scan)
#     t = requests.get(url).text
#     return parse_peaks(t)
#
# def parse_peaks(peak_string):
#     return [[float(p.split(" ")[0]),float(p.split(" ")[-1])] for p in peak_string.split("\n")[8:] if p != ""]
#
# def retreive_synthetics_w_peaks(protein,peptide):
#     query = {"overlaps":"false","libraries":"3,4","protein_name":protein,"sequence_input":peptide}
#     url = "https://massive.ucsd.edu/ProteoSAFe/ProteinCoverageServlet?task=protein_explorer_representatives&file=&pageSize=30&offset=0"
#     url += "&query={}".format(urllib.parse.quote(urllib.parse.quote('#' + json.dumps(query))))
#     url += "&query_type=representative&_=1558482309060"
#     results = requests.get(url).json()['row_data']
#     for result in results:
#         result['peaks'] = get_peaks(result['filename'],result['scan'])
#     return results
# would the psm be above the score threshold

def window_filter_peaks(peaks, window_size, top_peaks):
    peak_list_window_map = defaultdict(list)
    for peak in peaks:
        mass = peak[0]
        mass_bucket = int(mass/window_size)
        peak_list_window_map[mass_bucket].append(peak)

    new_peaks = []
    for bucket in peak_list_window_map:
        peaks_sorted_by_intensity = sorted(peak_list_window_map[bucket], key=lambda peak: peak[1], reverse=True)
        peaks_to_keep = peaks_sorted_by_intensity[:top_peaks]
        new_peaks += peaks_to_keep

    new_peaks = sorted(new_peaks, key=lambda peak: peak[0])
    return new_peaks

def precursor_filter_peaks(peaks, m, charge):
    new_peaks = []
    for peak in peaks:
        filtered = False
        for z in range(1,int(charge)+1):
            if abs(peak[0] - float(m)/z) < 2 or abs(peak[0] - (float(m)-18)/z) < 2:
                filtered = True
        if not filtered:
            new_peaks.append(peak)
    return new_peaks

def arguments():
    parser = argparse.ArgumentParser(description='mzTab to list of peptides')
    parser.add_argument('-m','--mztab', type = str, help='mzTab')
    parser.add_argument('-s','--spec_on_server', type = str, help='Spectra')
    parser.add_argument('-n','--novel_proteins', type = str, help='Novel Proteins')
    parser.add_argument('-p','--novel_psms', type = str, help='Novel PSMs')
    parser.add_argument('-e','--novel_peptides', type = str, help='Novel Peptides')
    parser.add_argument('-y','--synthetics', type = str, help='Matched Synthetics')
    parser.add_argument('-t','--peak_tolerance', type = str, help='Peak Tolerance for Matching')

    if len(sys.argv) < 4:
        parser.print_help()
        sys.exit(1)
    return parser.parse_args()

def main():
    args = arguments()
    ids = {}

    with open(args.synthetics,'rb') as f:
        synthetic_scans = pickle.load(f)

    tol = float(args.peak_tolerance)

    if args.mztab != "" and args.mztab != None:
        for mztab_file in glob.glob(args.mztab + '/*'):
            try:
                ids = mztab.read(mztab_file, ids)
            except:
                ids = mztab.read_lib(mztab_file, ids)
    peptide_to_psm = defaultdict(list)
    for filescan,psm in ids.items():
        peptide_to_psm[''.join([p.replace('I','L') for p in psm[0].sequence if p.isalpha()])].append((filescan,psm))

    peptide_to_protein = {}

    psms_to_consider = defaultdict(set)

    cosine_to_synthetic = defaultdict(lambda: (0,('N/A','N/A')))

    with open(args.novel_proteins) as f:
        r = csv.DictReader(f, delimiter = '\t')
        for l in r:
            supporting_peptides = []
            novel_peptides = []
            if l['supporting_peptides'] != ' ':
                supporting_peptides = [p.split(' ')[0] for p in l['supporting_peptides'].split(',')]
            if l['novel_peptides'] != ' ':
                novel_peptides = [p.split(' ')[0] for p in l['novel_peptides'].split(',')]
            protein = l['protein']
            for peptide in supporting_peptides + novel_peptides:
                for (filescan,psm) in peptide_to_psm[peptide.replace('I','L')]:
                    peptide_to_protein[peptide.replace('I','L')] = protein
                    if synthetic_scans.get((psm[0].sequence,psm[0].charge)):
                        psms_to_consider[filescan[0]].add(filescan[1])



    with open(args.novel_psms, 'w') as fw_psm, open(args.novel_peptides, 'w') as fw_pep:
        header = ['protein','filename','scan','sequence','charge','score','pass','type','parent_mass','synthetic_filename','synthetic_scan','best_cosine_to_synthetic','frag_tol']
        w_psm = csv.DictWriter(fw_psm, delimiter = '\t', fieldnames = header)
        w_pep = csv.DictWriter(fw_pep, delimiter = '\t', fieldnames = header)
        w_psm.writeheader()
        w_pep.writeheader()
        for peptide, protein in peptide_to_protein.items():
            best_psm = None
            best_psm_score = 0
            for (filescan,psm) in peptide_to_psm[peptide.replace('I','L')]:
                cosine, best_synthetic = cosine_to_synthetic[filescan]
                if best_synthetic[0] != 'N/A':
                    synthetic_filename = 'f.' + best_synthetic[0]
                else:
                    synthetic_filename = best_synthetic[0]
                psm_row = {
                    'protein':protein,
                    'filename':filescan[0],
                    'scan':filescan[1],
                    'sequence':psm[0].sequence,
                    'charge':psm[0].charge,
                    'score':psm[0].score,
                    'pass':'Above' if psm[0].score > CUTOFF[str(len(peptide))] else 'Below',
                    'parent_mass':psm[0].parent_mass,
                    'type':'Supporting' if peptide in supporting_peptides else 'Novel',
                    'synthetic_filename': synthetic_filename,
                    'synthetic_scan': best_synthetic[1],
                    'best_cosine_to_synthetic': cosine,
                    'frag_tol':psm[0].tolerance if psm[0].tolerance else 0
                }
                w_psm.writerow(psm_row)
                if psm[0].score > best_psm_score:
                    best_psm = psm_row
            w_pep.writerow(best_psm)

if __name__ == '__main__':
    main()
