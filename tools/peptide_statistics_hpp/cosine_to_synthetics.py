import argparse
import sys
import csv
from collections import defaultdict
from pyteomics import mzml, mzxml
from pathlib import Path
import pickle
from datetime import datetime
from python_ms_utilities import processing

CUTOFF = {"6": 2.12167006027, "7": 2.32154253572, "8": 2.53031824698, "9": 2.7717776307, "10": 3.0419577266, "11": 3.20404758311, "12": 3.4028353721, "13": 3.14854063781, "14": 3.47053358426, "15": 3.35416414915, "16": 3.34418455532, "17": 3.24904542906, "18": 3.18426919317, "19": 3.01362943461, "20": 3.12316407632, "21": 3.08160158158, "22": 2.59466460406, "23": 2.96230440256, "24": 3.11610830789, "25": 2.93990420679, "26": 2.6947192629, "27": 2.43042157531, "28": 2.5287628139300002, "29": 2.26034401643, "30": 2.60979068254, "31": 2.70004841417, "32": 2.69919517925, "33": 2.18110553715, "34": 1.90115111418, "35": 1.5402648112000001, "36": 1.74919562203, "37": 1.88066887473, "38": 1.58471929702, "39": 1.73377627878, "40": 3.10312899149}

def make_usi(filename, scan, sequence, charge):
    filename_path = Path(filename.replace('f.',''))
    dataset = filename_path.parts[0]
    file = filename_path.stem
    if 'MSV' in dataset or 'PXD' in dataset:
        usi = 'mzspec:{}:{}:scan:{}:{}/{}'.format(dataset, file, scan, sequence, charge)
    else:
        usi = 'N/A'

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

def get_mzml_spectrum(mzml_object, scan):
    #thermo only for now
    try:
        spectrum = mzml_object.get_by_id("controllerType=0 controllerNumber=1 scan={}".format(scan))
    except:
        spectrum = mzml_object.get_by_id("scan={}".format(scan))
    peaks = [processing.Peak(*p) for p in zip(spectrum['m/z array'],spectrum['intensity array'])]
    precursor = float(spectrum['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['selected ion m/z'])
    return processing.Spectrum(
        peaks,
        precursor,
        None,
        spectrum['id'],
        None
    )

def get_mzxml_spectrum(mzxml_object, scan):
    #thermo only for now
    spectrum = mzxml_object.get_by_id(scan)
    peaks = [processing.Peak(*p) for p in zip(spectrum['m/z array'],spectrum['intensity array'])]
    precursor = float(spectrum["precursorMz"][0]["precursorMz"])
    return processing.Spectrum(
        peaks,
        precursor,
        None,
        spectrum['id'],
        None
    )

def extract_annotated_peaks(spectrum, fragment_tolerance):
    spectrum = processing.filter_all_isobaric_tag_peaks(spectrum,fragment_tolerance)
    spectrum = processing.filter_precursor_peaks(spectrum,fragment_tolerance)
    spectrum = processing.window_filter_peaks(spectrum,50,10)
    spectrum = processing.normalize_spectrum(spectrum)
    explained_intensity, ion_vector = processing.calculate_explained_intensity(spectrum,fragment_tolerance)
    return explained_intensity, spectrum._replace(peaks = ion_vector)

def process_spectrum(psms_to_consider, filename, synthetic_scans, tol, threshold, peaks_obj, spectrum_select_func):
    cosine_to_synthetic = defaultdict(lambda: (-1,('N/A','N/A')))
    explained_intensity_per_spectrum = {}
    for scan in psms_to_consider[filename].keys():
        sequence = psms_to_consider[filename][scan]['sequence']
        charge = psms_to_consider[filename][scan]['charge']
        tolerance = psms_to_consider[filename][scan]['tolerance'] if psms_to_consider[filename][scan]['tolerance'] else tol
        matching_synthetics = synthetic_scans.get((sequence.replace('+229.163','').replace('+229.162932',''),charge),[])
        spectrum = spectrum_select_func(peaks_obj,scan)._replace(precursor_z = int(charge), annotation = processing.Annotation(sequence, None))
        spectrum_ei, spectrum_ion_vector = extract_annotated_peaks(spectrum, tolerance)
        explained_intensity_per_spectrum[(filename,scan)] = spectrum_ei
        if threshold == 0:
            for synthetic_filescan, _ in matching_synthetics:
                cosine_to_synthetic[(filename,scan)] = (0,synthetic_filescan)
        else:
            cosine = 0
            for synthetic_filescan,synthetic_ion_vector in matching_synthetics:
                cosine = processing.match_peaks(spectrum_ion_vector, synthetic_ion_vector, tolerance)
                if cosine > cosine_to_synthetic[(filename,scan)][0]:
                    cosine_to_synthetic[(filename,scan)] = (cosine,synthetic_filescan)
    return cosine_to_synthetic, explained_intensity_per_spectrum

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
            synthetic_keys.add((l['sequence'].replace('+229.163','').replace('+229.162932',''),l['charge']))
            all_psms.append(l)
            try:
                tolerance = float(l['frag_tol'])
            except:
                tolerance = None
            psms_to_consider[l['filename']][l['scan']] = {'sequence':l['sequence'],'charge':l['charge'], 'tolerance':tolerance}
    print("{}: Finished reading psms".format(datetime.now().strftime("%H:%M:%S")))


    tol = float(args.peak_tolerance)
    threshold = float(args.cosine_threshold)
    explained_intensity = float(args.explained_intensity)

    if len(synthetic_keys) > 0:
        print("{}: Loading synthetics".format(datetime.now().strftime("%H:%M:%S")))
        with open(args.synthetics,'rb') as f:
            total_subbuckets = pickle.load(f)
            print("Total buckets: {}:".format(total_subbuckets))
            for _ in range(total_subbuckets):
                partial_loaded_synthetics = pickle.load(f)
                print("Bucket size: {}:".format(len(partial_loaded_synthetics.keys())))
                for key in synthetic_keys:
                    if key in partial_loaded_synthetics:
                        for synthetic_filescan,synthetic_spectrum in partial_loaded_synthetics[key]:
                            _, synthetic_ion_vector = extract_annotated_peaks(
                                synthetic_spectrum._replace(
                                    precursor_z = int(key[1]),
                                    annotation = processing.Annotation(key[0], None)
                                ), tolerance
                            )
                            synthetic_scans[key].append((synthetic_filescan,synthetic_ion_vector))
                print("Cumulative total: {}:".format(len(synthetic_scans.keys())))
        print("{}: Loaded {} synthetics".format(datetime.now().strftime("%H:%M:%S"),len(synthetic_scans)))
    else:
        print("Not loading synthetics, nothing to match")

    cosine_to_synthetic = defaultdict(lambda: (-1,('N/A','N/A')))
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

        if threshold > 0 or explained_intensity > 0:
            if not Path(filepath).exists():
                print("File {} likely moved or doesn't exist.".format(filepath))
            else:
                print("{}: About to read {}".format(datetime.now().strftime("%H:%M:%S"),filename))

                file = None
                file_object = None
                get_spectrum_func = None

                if exts[0] == '.mzML':
                    file = open(filepath, 'rb')
                    file_object = mzml.PreIndexedMzML(file)
                    get_spectrum_func = get_mzml_spectrum
                if exts[0] == '.mgf':
                    pass
                if exts[0] == '.mzXML':
                    file = open(filepath, 'rb')
                    file_object = mzxml.MzXML(file,use_index = True)
                    get_spectrum_func = get_mzxml_spectrum

                if file:
                    if file_object and get_spectrum_func:
                        file_cosine_to_synthetic, file_explained_intensity_per_spectrum = process_spectrum(psms_to_consider,filename,synthetic_scans,tol,threshold,file_object,get_spectrum_func)
                        cosine_to_synthetic.update(file_cosine_to_synthetic)
                        explained_intensity_per_spectrum.update(file_explained_intensity_per_spectrum)
                    file.close()
        else:
            for scan in psms_to_consider[filename].keys():
                sequence = psms_to_consider[filename][scan]['sequence']
                charge = psms_to_consider[filename][scan]['charge']
                matching_synthetics = synthetic_scans.get((sequence.replace('+229.163','').replace('+229.162932',''),charge),[])
                for synthetic_filescan, _ in matching_synthetics:
                    cosine_to_synthetic[(filename,scan)] = (0,synthetic_filescan)

    print("{}: About to write out PSMs".format(datetime.now().strftime("%H:%M:%S")))

    with open(args.output_psms.joinpath(args.jobs.name), 'w') as fw_psm:
        header = psms_header + ['usi','synthetic_filename','synthetic_scan','synthetic_usi','cosine','explained_intensity']
        w_psm = csv.DictWriter(fw_psm, delimiter = '\t', fieldnames = header)
        w_psm.writeheader()
        for psm in all_psms:
            cosine, best_synthetic = cosine_to_synthetic[(psm['filename'],psm['scan'])]
            if best_synthetic[0] != 'N/A':
                synthetic_filename = 'f.' + best_synthetic[0].replace('/data/massive/','')
                synthetic_scan = best_synthetic[1]
            else:
                synthetic_filename = best_synthetic[0]
                synthetic_scan = best_synthetic[1]
            psm['usi'] = make_usi(psm['filename'], psm['scan'], psm['sequence'], psm['charge'])
            psm['synthetic_filename'] = synthetic_filename
            psm['synthetic_scan'] = synthetic_scan
            psm['synthetic_usi'] = make_usi(synthetic_filename, synthetic_scan, psm['sequence'].replace('+229.163','').replace('+229.162932',''), psm['charge'])
            psm['cosine'] = cosine
            psm['explained_intensity'] = explained_intensity_per_spectrum.get((psm['filename'],psm['scan']),0)
            psm['frag_tol'] = psm['frag_tol'] if psm['frag_tol'] != 0 else tol
            w_psm.writerow(psm)
    print("{}: Finished writing out PSMs".format(datetime.now().strftime("%H:%M:%S")))

if __name__ == '__main__':
    main()
