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


def create_peaks_url(filename,scan):
    filename = urllib.parse.quote(">" + filename)
    url = "http://ccms-internal.ucsd.edu/ProteoSAFe/DownloadResultFile"
    url += "?invoke=annotatedSpectrumImageText"
    url += "&task=protein_explorer_provenance"
    url += "&block=0"
    url += "&file=FILE-" + filename
    url += "&peptide=*..*"
    url += "&scan=" + str(scan)
    url += "&uploadfile=True&task=4f2ac74ea114401787a7e96e143bb4a1"
    return url

def get_peaks(filename, scan):
    url = create_peaks_url(filename,scan)
    t = requests.get(url).text
    return parse_peaks(t)

def parse_peaks(peak_string):
    return [[float(p.split(" ")[0]),float(p.split(" ")[-1])] for p in peak_string.split("\n")[8:] if p != ""]

def retreive_synthetics_w_peaks(protein,peptide):
    query = {"overlaps":"false","libraries":"3,4","protein_name":protein,"sequence_input":peptide}
    url = "https://massive.ucsd.edu/ProteoSAFe/ProteinCoverageServlet?task=protein_explorer_provenance&file=&pageSize=30&offset=0"
    url += "&query={}".format(urllib.parse.quote(urllib.parse.quote('#' + json.dumps(query))))
    url += "&query_type=provenance&_=1558482309060"
    results = requests.get(url).json()['row_data']
    for result in results:
        result['peaks'] = get_peaks(result['filename'],result['scan'])
    return results
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
        header = ['protein','filename','scan','sequence','charge','type','parent_mass','synthetic_filename','synthetic_scan','best_cosine_to_synthetic']
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
                synthetic_psms = retreive_synthetics_w_peaks(protein, peptide)
                max_cosine = 0
                best_synthetic = ('N/A','N/A')
                for (filescan,psm) in peptide_to_psm[peptide]:
                    peaks = get_peaks(filescan[0], filescan[1])
                    for synthetic in synthetic_psms:
                        if int(synthetic['charge']) == psm[0].charge:
                            score,_ = sa.score_alignment(peaks,synthetic['peaks'],0,0,0.05)
                            if (score > max_cosine):
                                max_cosine = score
                                best_synthetic = (synthetic['filename'],synthetic['scan'])
                    w.writerow({
                        'protein':protein,
                        'filename':filescan[0],
                        'scan':filescan[1],
                        'sequence':psm[0].sequence,
                        'charge':psm[0].charge,
                        'parent_mass':psm[0].parent_mass,
                        'type':'Supporting' if peptide in supporting_peptides else 'Novel',
                        'synthetic_filename': best_synthetic[0],
                        'synthetic_scan': best_synthetic[1],
                        'best_cosine_to_synthetic': max_cosine
                    })



if __name__ == '__main__':
    main()
