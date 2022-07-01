import argparse
import sys
import csv
from collections import defaultdict
from pyteomics import mzml, mzxml, mgf
from pathlib import Path
import pickle
from datetime import datetime
from python_ms_utilities import processing

def PathList(argument_string):
    output_pathlist = []
    if argument_string is None:
        pass
    elif argument_string.rstrip() == '':
        pass
    elif Path(argument_string).exists():
        input_path = Path(argument_string)
        if input_path.is_dir():
            output_pathlist = list(input_path.glob('*'))
        else:
            output_pathlist = [input_path]
    return output_pathlist

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
    parser.add_argument('-p','--spectrum_files', type = Path, help='Spectrum Files (if not public)')
    parser.add_argument('-t','--peak_tolerance', type = str, help='Peak Tolerance for Matching')
    parser.add_argument('-e','--explained_intensity', type = str, help='Explained Intensity Filter')
    parser.add_argument('-l','--labeled', type = str, help='Labeled data')
    parser.add_argument('-c','--cosine_threshold', type = str, help='Cosine Threshold')
    parser.add_argument('-m','--low_mass_filter', type = float, help='Minimum m/z to retain', default=232)
    parser.add_argument('-s','--min_snr', type = float, help='SNR threshold', default=2)
    parser.add_argument('-f','--filter_synthetics', type = str, help='Filter Matched Synthetics', default='No')

    if len(sys.argv) < 4:
        parser.print_help()
        sys.exit(1)
    return parser.parse_args()

def read_mzml_spectrum(spectrum):
    peaks = [processing.Peak(float(p[0]),float(p[1])) for p in zip(spectrum['m/z array'],spectrum['intensity array'])]
    precursor = float(spectrum['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['selected ion m/z'])
    return processing.Spectrum(
        peaks,
        precursor,
        None,
        spectrum['id'],
        None
    )

def get_mzml_spectrum(mzml_object, scan):
    #thermo only for now
    try:
        spectrum = mzml_object.get_by_id("controllerType=0 controllerNumber=1 scan={}".format(scan))
    except:
        spectrum = mzml_object.get_by_id("scan={}".format(scan))
    return read_mzml_spectrum(spectrum)

def read_mzxml_spectrum(spectrum):
    peaks = [processing.Peak(float(p[0]),float(p[1])) for p in zip(spectrum['m/z array'],spectrum['intensity array'])]
    precursor = float(spectrum["precursorMz"][0]["precursorMz"])
    return processing.Spectrum(
        peaks,
        precursor,
        None,
        "scan={}".format(spectrum['id']),
        None
    )

def get_mzxml_spectrum(mzxml_object, scan):
    #thermo only for now
    spectrum = mzxml_object.get_by_id(scan)
    return read_mzxml_spectrum(spectrum)

def extract_annotated_peaks(spectrum, fragment_tolerance, low_mass_filter, min_snr):
    spectrum, explained_intensity, ion_vector, b_y_peaks = processing.process_spectrum(
            spectrum,
            fragment_tolerance,
            precursor_filter_window=1.5,
            low_mass_filter=low_mass_filter,
            isobaric_tag_type=None,
            min_snr=min_snr,
            num_top_unannotated_envelopes_to_remove=2
    )
    ion_vector = spectrum._replace(peaks = ion_vector)
    ion_vector = processing.normalize_spectrum(ion_vector)
    return (explained_intensity,b_y_peaks), ion_vector

def find_ei_and_intensity(spectrum, psm, synthetic_scans, tol, low_mass_filter, min_snr):
    best_cosine = None
    sequence = psm['sequence']
    charge = psm['charge']
    tolerance = psm['tolerance'] if psm['tolerance'] else tol
    spectrum = spectrum._replace(precursor_z = int(charge), annotation = processing.Annotation(sequence, None))
    matching_synthetics = synthetic_scans.get((sequence.replace('+229.163','').replace('+229.162932',''),charge),[])
    spectrum_ei, spectrum_ion_vector = extract_annotated_peaks(spectrum, tolerance, low_mass_filter, min_snr)
    for synthetic_filescan, synthetic_ion_vector in matching_synthetics:
        cosine = processing.match_peaks(spectrum_ion_vector, synthetic_ion_vector, tolerance)
        if not best_cosine or cosine > best_cosine[0]:
            best_cosine = (cosine,synthetic_filescan)
    return spectrum_ei, best_cosine

def process_spectrum(psms_to_consider, filename, synthetic_scans, tol, low_mass_filter, min_snr, threshold, peaks_obj, spectrum_select_func):
    cosine_to_synthetic = defaultdict(lambda: (-1,('N/A','N/A')))
    explained_intensity_per_spectrum = {}
    for scan in psms_to_consider[filename].keys():
        spectrum = spectrum_select_func(peaks_obj,scan)
        spectrum_ei, cosine_w_file = find_ei_and_intensity(spectrum,psms_to_consider[filename][scan],synthetic_scans, tol, low_mass_filter, min_snr)
        explained_intensity_per_spectrum[(filename,scan)] = spectrum_ei
        if cosine_w_file:
            cosine_to_synthetic[(filename,scan)] = cosine_w_file
    return cosine_to_synthetic, explained_intensity_per_spectrum

def process_spectrum_read_file(psms_to_consider, filename, synthetic_scans, tol, low_mass_filter, min_snr, threshold, reader, spectrum_select_func,read_scan):
    start_time = datetime.now()
    cosine_to_synthetic = defaultdict(lambda: (-1,('N/A','N/A')))
    explained_intensity_per_spectrum = {}
    all_spectra = []
    for i,s in enumerate(reader):
        if (i != 0 and i%10000 == 0):
            elapsed_time = (datetime.now()-start_time).seconds
            rate = int(i/elapsed_time) if elapsed_time != 0 else int(i)
            print("{}: Read {} scans ({}/second)".format(datetime.now().strftime("%H:%M:%S"),i,rate))
        all_spectra.append(s)

    start_time = datetime.now()
    spectra_to_process = 0
    for s in all_spectra:
        scan = read_scan(s)
        if scan in psms_to_consider[filename]:
            if (spectra_to_process != 0 and spectra_to_process%100 == 0):
                elapsed_time = (datetime.now()-start_time).seconds
                rate = int(spectra_to_process/elapsed_time) if elapsed_time != 0 else int(spectra_to_process)
                print("{}: Processed {} scans ({}/second)".format(datetime.now().strftime("%H:%M:%S"),spectra_to_process,rate))
            spectrum = spectrum_select_func(s)
            spectrum_ei, cosine_w_file = find_ei_and_intensity(spectrum,psms_to_consider[filename][scan],synthetic_scans, tol, low_mass_filter, min_snr)
            explained_intensity_per_spectrum[(filename,scan)] = spectrum_ei
            if cosine_w_file:
                cosine_to_synthetic[(filename,scan)] = cosine_w_file
            spectra_to_process += 1

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

    if len(synthetic_keys) > 0 and args.synthetics:
        start_synthetic_read = datetime.now()
        print("{}: Loading synthetics".format(start_synthetic_read.strftime("%H:%M:%S")))
        synthetics_loaded = 0
        for synthetics_path in args.synthetics.glob('*'):
            with open(synthetics_path) as synthetics_file:
                with mgf.read(synthetics_file) as reader:
                    for i,s in enumerate(reader):
                        peptide = s['params'].get('seq')
                        charge = int(s['params'].get('charge')[0])
                        # print(peptide, str(charge))
                        if (peptide, str(charge)) in synthetic_keys:
                            synthetics_loaded += 1
                            filename = s['params'].get('originalfile_filename',s['params'].get('provenance_filename'))
                            scan = s['params'].get('originalfile_scan',s['params'].get('provenance_scan'))
                            precursor_mz = s['params'].get('pepmass')[0]
                            mz = s['m/z array']
                            intn = s['intensity array']
                            peaks = [processing.Peak(float(mz[i]),float(intn[i])) for i in range(len(mz))]
                            spectrum = processing.Spectrum(
                                peaks,
                                precursor_mz,
                                charge,
                                None,
                                processing.Annotation(peptide, None)
                            )
                            synthetic_low_mass, synthetic_min_snr = (args.low_mass_filter, args.min_snr) if args.filter_synthetics == 'Yes' else (0,0)
                            _, synthetic_ion_vector = extract_annotated_peaks(spectrum, 0.05, synthetic_low_mass, synthetic_min_snr)
                            synthetic_scans[(peptide, str(charge))].append(((filename,scan),synthetic_ion_vector))
        print("{}: Loaded {} synthetics".format(datetime.now().strftime("%H:%M:%S"),synthetics_loaded))
    else:
        print("Not loading synthetics")

    cosine_to_synthetic = defaultdict(lambda: (-1,('N/A','N/A')))
    explained_intensity_per_spectrum = {}

    min_spectra_to_load_file = 20

    for filename in psms_to_consider:
        print("{}: Looking at {}".format(datetime.now().strftime("%H:%M:%S"),filename))
        exts = Path(filename).suffixes
        if 'MSV' in filename:
            filepath = filename.replace('f.','/data/massive/')
            if filepath[0] == 'M':
                filepath = '/data/massive/' + filepath
        elif 'f.' in filename:
            # checking uploads directory
            filepath = filename.replace('f.','/data/ccms-data/uploads/')
        else:
            filepath = filename

        file_exists = False


        potential_matches = []
        if not Path(filepath).exists():
            for potential_file in args.spectrum_files.glob('**'):
                if potential_file.name == Path(filename).name:
                    potential_matches.append(potential_file)
        
        if len(potential_matches) > 1:
            depth = -1
            filename_parts = Path(filename).parts
            while len(potential_matches) > 1:
                potential_match_update = []
                for match in potential_matches:
                    if match.parts[depth] == filename_parts[depth]:
                        potential_match_update.append(match)
                depth -= 1
                potential_matches = potential_match_update
            filepath = potential_matches[0]
        elif len(potential_matches) == 1:
            filepath = potential_matches[0]

        if not Path(filepath).exists():
            print("File {} likely moved or doesn't exist, please upload the file in the input file.".format(filepath))
        else:
            print("{}: About to read {} ({} PSMs)".format(datetime.now().strftime("%H:%M:%S"),filename,len(psms_to_consider[filename])))
            file_exists = True

        if threshold > 0 or explained_intensity > 0:
            if file_exists:

                file = None
                file_object = None
                get_spectrum_func = None
                read_scan = None

                if len(psms_to_consider[filename]) >= min_spectra_to_load_file:
                    if exts[0] == '.mzML':
                        file = open(filepath, 'rb')
                        file_object = mzml.MzML(file)
                        get_spectrum_func = read_mzml_spectrum
                        read_scan = lambda s: str(s['id'].split('scan=')[1])
                    if exts[0] == '.mgf':
                        pass
                    if exts[0] == '.mzXML':
                        file = open(filepath, 'rb')
                        file_object = mzxml.MzXML(file)
                        get_spectrum_func = read_mzxml_spectrum
                        read_scan = lambda s: str(s['id'])
                else:
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
                        print("{}: Opened file {}, ready to load".format(datetime.now().strftime("%H:%M:%S"),filename))
                        if len(psms_to_consider[filename]) >= min_spectra_to_load_file:
                            file_cosine_to_synthetic, file_explained_intensity_per_spectrum = process_spectrum_read_file(psms_to_consider,filename,synthetic_scans,tol,args.low_mass_filter,args.min_snr,threshold,file_object,get_spectrum_func,read_scan)
                        else:
                            file_cosine_to_synthetic, file_explained_intensity_per_spectrum = process_spectrum(psms_to_consider,filename,synthetic_scans,tol,args.low_mass_filter,args.min_snr,threshold,file_object,get_spectrum_func)
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
        header = psms_header + ['usi','synthetic_filename','synthetic_scan','synthetic_usi','cosine','explained_intensity','matched_ions']
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
            ei, num_matched_peaks = explained_intensity_per_spectrum.get((psm['filename'],psm['scan']),(0,0))
            psm['explained_intensity'] = ei
            psm['matched_ions'] = num_matched_peaks

            psm['frag_tol'] = psm['frag_tol'] if psm['frag_tol'] != 0 else tol
            w_psm.writerow(psm)
    print("{}: Finished writing out PSMs".format(datetime.now().strftime("%H:%M:%S")))

if __name__ == '__main__':
    main()
