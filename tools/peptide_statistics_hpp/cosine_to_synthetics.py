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
from pyteomics import mzml, mzxml
import ming_spectrum_library as msl
import ming_psm_library as mpl
from pathlib import Path
import pickle
from datetime import datetime




CUTOFF = {"6": 2.12167006027, "7": 2.32154253572, "8": 2.53031824698, "9": 2.7717776307, "10": 3.0419577266, "11": 3.20404758311, "12": 3.4028353721, "13": 3.14854063781, "14": 3.47053358426, "15": 3.35416414915, "16": 3.34418455532, "17": 3.24904542906, "18": 3.18426919317, "19": 3.01362943461, "20": 3.12316407632, "21": 3.08160158158, "22": 2.59466460406, "23": 2.96230440256, "24": 3.11610830789, "25": 2.93990420679, "26": 2.6947192629, "27": 2.43042157531, "28": 2.5287628139300002, "29": 2.26034401643, "30": 2.60979068254, "31": 2.70004841417, "32": 2.69919517925, "33": 2.18110553715, "34": 1.90115111418, "35": 1.5402648112000001, "36": 1.74919562203, "37": 1.88066887473, "38": 1.58471929702, "39": 1.73377627878, "40": 3.10312899149}

def make_usi(filename, scan, sequence, charge):
    filename_path = Path(filename.replace('f.',''))
    dataset = filename_path.parts[0]
    file = filename_path.stem
    if 'MSV' in dataset or 'PXD' in dataset:
        usi = 'mzspec:{}:{}:scan:{}:{}/{}'.format(dataset, file, scan, sequence, charge)
    else:
        usi = None

    return usi

def arguments():
    parser = argparse.ArgumentParser(description='mzTab to list of peptides')
    parser.add_argument('-y','--synthetics', type = Path, help='Matched Synthetics')
    parser.add_argument('-j','--jobs', type = Path, help='Jobs for parallelism')
    parser.add_argument('-i','--input_psms', type = Path, help='Input PSMs')
    parser.add_argument('-o','--output_psms', type = Path, help='Output PSMs')
    parser.add_argument('-t','--peak_tolerance', type = str, help='Peak Tolerance for Matching')
    parser.add_argument('-e','--explained_intensity', type = str, help='Explained Intensity Filter')
    parser.add_argument('-l','--labeled', type = str, help='Labeled data')

    parser.add_argument('-c','--cosine_threshold', type = str, help='Cosine Threshold')

    if len(sys.argv) < 4:
        parser.print_help()
        sys.exit(1)
    return parser.parse_args()

def prepare_spectrum(peaks, tol, precursor, charge, sequence, process = False):

    #filter TMT
    msl.filter_isobaric_label_peaks(peaks,tol)

    charge_set = range(1, int(charge))

    theoretical_peaks = mpl.create_theoretical_peak_map(sequence, ["b",  "b-iso", "y", "y-iso", "b-H2O", "b-NH3", "y-H2O", "y-NH3", "a"], charge_set=charge_set)

    if process:
        peaks = msl.filter_precursor_peaks(peaks,2,precursor)
        peaks = msl.window_filter_peaks(peaks, 50, 10)

    annotated_peaks, unannotated_peaks, ion_vector = mpl.extract_annotated_peaks(theoretical_peaks, peaks, tol)

    peaks = sa.normalize_spectrum(sa.convert_to_peaks(peaks))

    explained_intensity = sum(p[0] for p in annotated_peaks)/sum(p[0] for p in annotated_peaks+unannotated_peaks)

    annotated_peaks = sa.normalize_spectrum(sa.convert_to_peaks(annotated_peaks))
    ion_vector = sa.normalize_spectrum(sa.convert_to_peaks(ion_vector))
    # print("Peaks ", peaks)
    # print("Annotated Peaks ", annotated_peaks)

    return peaks, annotated_peaks, ion_vector, explained_intensity

def main():
    args = arguments()

    synthetic_scans = defaultdict(list)

    synthetic_keys = set()

    psms_to_consider = defaultdict(lambda: defaultdict(dict))
    all_psms = []
    psms_header = []

    print("{}: Reading psms".format(datetime.now().strftime("%H:%M:%S")))
    with open(args.jobs) as f:
        r = csv.DictReader(f, delimiter='\t')
        psms_header = r.fieldnames
        for l in r:
            synthetic_keys.add((l['sequence'].replace('+229.163',''),l['charge']))
            all_psms.append(l)
            psms_to_consider[l['filename']][l['scan']] = {'sequence':l['sequence'].replace('+229.163',''),'charge':l['charge']}
    print("{}: Finished reading psms".format(datetime.now().strftime("%H:%M:%S")))

    if len(synthetic_keys) > 0:
        print("{}: Loading synthetics".format(datetime.now().strftime("%H:%M:%S")))
        with open(args.synthetics,'rb') as f:
            for _ in range(pickle.load(f)):
                partial_loaded_synthetics = pickle.load(f)
                for key in synthetic_keys:
                    if key in partial_loaded_synthetics:
                        synthetic_scans[key].extend(partial_loaded_synthetics[key])
        print("{}: Loaded {} synthetics".format(datetime.now().strftime("%H:%M:%S"),len(synthetic_scans)))
    else:
        print("Not loading synthetics, nothing to match")

    tol = float(args.peak_tolerance)
    threshold = float(args.cosine_threshold)

    cosine_to_synthetic = defaultdict(lambda: (0,('N/A','N/A')))
    cosine_to_synthetic_annotated = defaultdict(lambda: (0,('N/A','N/A')))
    explained_intensity_per_spectrum = {}

    for filename in psms_to_consider:
        print("{}: Looking at {}".format(datetime.now().strftime("%H:%M:%S"),filename))
        exts = Path(filename).suffixes
        if 'MSV' in filename:
            filepath = filename.replace('f.','/data/massive/')
            if filepath[0] == 'M':
                filepath = '/data/massive/' + filepath
        else:
            filepath = filename
        if not Path(filepath).exists():
            print("File {} likely moved or doesn't exist.".format(filepath))
        else:
            if exts[0] == '.mzML':
                precursor_func = lambda spectrum: float(spectrum['precursorList']['precursor'][0]['isolationWindow']['isolation window target m/z'])
                with open(filepath, 'rb') as f:
                    mzml_object = mzml.PreIndexedMzML(f)
                    for scan in psms_to_consider[filename].keys():
                        sequence = psms_to_consider[filename][scan]['sequence']
                        charge = psms_to_consider[filename][scan]['charge']
                        matching_synthetics = synthetic_scans.get((sequence.replace('+229.163',''),charge),[])
                        spectrum = mzml_object.get_by_id("controllerType=0 controllerNumber=1 scan={}".format(scan))
                        peaks = list(zip(spectrum['m/z array'],spectrum['intensity array']))
                        peaks, annotated_peaks, ion_vector, explained_intensity = prepare_spectrum(peaks,tol,precursor_func(spectrum),charge,sequence, True)
                        explained_intensity_per_spectrum[(filename,scan)] = explained_intensity
                        if threshold == 0:
                            for synthetic_filescan, synthetic_peaks in matching_synthetics:
                                cosine_to_synthetic[(filename,scan)] = (0,synthetic_filescan)
                                cosine_to_synthetic_annotated[(filename,scan)] = (0,synthetic_filescan)

                        else:
                            # print("{}: About to calculate {} cosines".format(datetime.now().strftime("%H:%M:%S"),len(matching_synthetics)))
                            for synthetic_filescan,synthetic_peaks in matching_synthetics:
                                synthetic_peaks, annotated_synthetic_peaks, ion_vector_synthetic, _ = prepare_spectrum(synthetic_peaks,tol,precursor_func(spectrum),charge,sequence.replace('+229.163',''))
                                cosine,_ = sa.score_alignment(peaks,synthetic_peaks,0,0,tol)
                                cosine_annotated,_ = sa.score_alignment(ion_vector,ion_vector_synthetic,0,0,tol)
                                if cosine > cosine_to_synthetic[(filename,scan)][0]:
                                    cosine_to_synthetic[(filename,scan)] = (cosine,synthetic_filescan)
                                if cosine_annotated > cosine_to_synthetic_annotated[(filename,scan)][0]:
                                    cosine_to_synthetic_annotated[(filename,scan)] = (cosine_annotated,synthetic_filescan)
                # print("{}: About to read {}".format(datetime.now().strftime("%H:%M:%S"),filename))
                # with mzml.read(filepath) as reader:
                #     precursor_func = lambda spectrum: float(spectrum['precursorList']['precursor'][0]['isolationWindow']['isolation window target m/z'])
                #     for spectrum in reader:
                #         spec_id = {idx.split('=')[0]:idx.split('=')[1] for idx in spectrum['id'].split(' ')}
                #         scan = spec_id['scan']
                #         if scan in psms_to_consider[filename]:
                #             peaks = zip(spectrum['m/z array'],spectrum['intensity array'])
                #             sequence = psms_to_consider[filename][scan]['sequence']
                #             charge = psms_to_consider[filename][scan]['charge']
                #             print((sequence,charge))
                #             matching_synthetics = synthetic_scans.get((sequence,charge),[])
                #             if len(matching_synthetics) > 0:
                #                 peaks = prepare_spectrum(peaks,tol,precursor_func(spectrum))
                #                 print("{}: {}".format(datetime.now().strftime("%H:%M:%S"),peaks))
                #             print("{}: About to calculate {} cosines".format(datetime.now().strftime("%H:%M:%S"),len(matching_synthetics)))
                #             for synthetic_filescan,synthetic_peaks in matching_synthetics:
                #                 cosine,_ = sa.score_alignment(peaks,synthetic_peaks,0,0,tol)
                #                 if cosine > cosine_to_synthetic[(filename,scan)][0]:
                #                     cosine_to_synthetic[(filename,scan)] = (cosine,synthetic_filescan)
            if exts[0] == '.mgf':
                pass
            if exts[0] == '.mzXML':
                print("{}: About to read {}".format(datetime.now().strftime("%H:%M:%S"),filename))
                precursor_func = lambda spectrum: float(spectrum["precursorMz"][0]["precursorMz"])
                with open(filepath, 'rb') as f:
                    mzxml_object = mzxml.MzXML(f,use_index = True)
                    # for spectrum in reader:
                    for scan in psms_to_consider[filename].keys():
                        sequence = psms_to_consider[filename][scan]['sequence']
                        charge = psms_to_consider[filename][scan]['charge']
                        matching_synthetics = synthetic_scans.get((sequence.replace('+229.163',''),charge),[])
                        spectrum = mzxml_object.get_by_id(scan)
                        peaks = list(zip(spectrum['m/z array'],spectrum['intensity array']))
                        peaks, annotated_peaks, ion_vector, explained_intensity = prepare_spectrum(peaks,tol,precursor_func(spectrum),charge,sequence, True)
                        explained_intensity_per_spectrum[(filename,scan)] = explained_intensity
                        if threshold == 0:
                            for synthetic_filescan, synthetic_peaks in matching_synthetics:
                                cosine_to_synthetic[(filename,scan)] = (0,synthetic_filescan)
                                cosine_to_synthetic_annotated[(filename,scan)] = (0,synthetic_filescan)
                        else:
                            # print("{}: About to calculate {} cosines".format(datetime.now().strftime("%H:%M:%S"),len(matching_synthetics)))
                            for synthetic_filescan,synthetic_peaks in matching_synthetics:
                                synthetic_peaks, annotated_synthetic_peaks, ion_vector, _ = prepare_spectrum(synthetic_peaks,tol,precursor_func(spectrum),charge,sequence.replace('+229.163',''))
                                cosine,_ = sa.score_alignment(peaks,synthetic_peaks,0,0,tol)
                                cosine_annotated,_ = sa.score_alignment(annotated_peaks,annotated_synthetic_peaks,0,0,tol)
                                if cosine > cosine_to_synthetic[(filename,scan)][0]:
                                    cosine_to_synthetic[(filename,scan)] = (cosine,synthetic_filescan)
                                if cosine_annotated > cosine_to_synthetic_annotated[(filename,scan)][0]:
                                    cosine_to_synthetic_annotated[(filename,scan)] = (cosine_annotated,synthetic_filescan)

    print("{}: About to write out PSMs".format(datetime.now().strftime("%H:%M:%S")))

    with open(args.output_psms.joinpath(args.jobs.name), 'w') as fw_psm:
        header = psms_header + ['usi','synthetic_filename','synthetic_scan','synthetic_usi','cosine','explained_intensity']
        w_psm = csv.DictWriter(fw_psm, delimiter = '\t', fieldnames = header)
        w_psm.writeheader()
        for psm in all_psms:
            cosine, best_synthetic = cosine_to_synthetic_annotated[(psm['filename'],psm['scan'])]
            if best_synthetic[0] != 'N/A':
                synthetic_filename = 'f.' + best_synthetic[0].replace('/data/massive/','')
                synthetic_scan = best_synthetic[1]
            else:
                synthetic_filename = best_synthetic[0]
                synthetic_scan = best_synthetic[1]
            psm['usi'] = make_usi(psm['filename'], psm['scan'], psm['sequence'], psm['charge'])
            psm['synthetic_filename'] = synthetic_filename
            psm['synthetic_scan'] = synthetic_scan
            psm['synthetic_usi'] = make_usi(synthetic_filename, synthetic_scan, psm['sequence'].replace('+229.163',''), psm['charge'])
            psm['cosine'] = cosine
            psm['explained_intensity'] = explained_intensity_per_spectrum.get((psm['filename'],psm['scan']),0)
            w_psm.writerow(psm)
    print("{}: Finished writing out PSMs".format(datetime.now().strftime("%H:%M:%S")))

if __name__ == '__main__':
    main()
