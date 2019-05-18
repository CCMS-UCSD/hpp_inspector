from collections import defaultdict
from csv import DictReader, DictWriter
from quantlib import mass

def read_ids(max_quant_ids, file_mapping):
    all_ids = []
    feature_to_msms = defaultdict(list)
    ids = defaultdict(list)
    ids_list = []
    error = []
    with open(max_quant_ids) as f:
        r = DictReader(f, delimiter = '\t')
        for l in r:
#             if l['Type'] != 'SECPEP':
            demangled_file = l['Raw file']
            if "spectrum-" in l['Raw file']:
                demangled_file = file_mapping[l['Raw file']].split('/')[-1].split('.')[0]
            pep_mz = calculate_theoretica_mz(format_mq_pep(l['Modified sequence']),l['Charge'])
            ms2_mz = float(l['m/z'])
            charge = int(l['Charge'])
            ppm_difference = min(list(enumerate([
                abs(ms2_mz - pep_mz - 4/charge),
                abs(ms2_mz - pep_mz - 3/charge),
                abs(ms2_mz - pep_mz - 2/charge),
                abs(ms2_mz - pep_mz - 1/charge),
                abs(ms2_mz - pep_mz),
                abs(ms2_mz - pep_mz + 1/charge),
                abs(ms2_mz - pep_mz + 2/charge),
                abs(ms2_mz - pep_mz + 3/charge),
                abs(ms2_mz - pep_mz + 4/charge)
            ])), key = lambda x: x[1])
            isotope = ppm_difference[0] - 4
#             isotopes.append(isotope)
            ppm_difference = (ppm_difference[1]/pep_mz)*1e6
            error.append(ppm_difference)
            if l['Type'] != 'MULTI-SECPEP':
                feature_to_msms[int(l['Evidence ID'])].append((demangled_file,int(l['Scan number']),format_mq_pep(l['Modified sequence']),l['Proteins'],l['Charge']))
                ids[(demangled_file,int(l['Scan number']))].append((format_mq_pep(l['Modified sequence']),l['Proteins'],l['Charge'],l['Type'],l['Score']))
                ids_list.append((demangled_file,int(l['Scan number']),format_mq_pep(l['Modified sequence']),l['Proteins'],l['Charge']))
    sns.distplot(error, kde = False)
    return ids, ids_list, feature_to_msms

def read_ms_scans(max_quant_scans):
    scan_to_rt = {}
    with open(max_quant_scans) as f:
        r = DictReader(f, delimiter = '\t')
        for l in r:
            scan_to_rt[(l['Raw file'],int(l['Scan number']))] = float(l['Retention time'])
    return scan_to_rt

def read_evidence(max_quant_evdience, file_mapping):
    features = []
    peptides = {}
    peptide_count = 0
    with open(max_quant_evdience) as f:
        r = DictReader(f, delimiter = '\t')
        for i,l in enumerate(r):
            l['Raw file'] = file_mapping[l['Raw file']]
            #modified_sequence = mass.format_mq_pep(l['Modified sequence'])
            #if modified_sequence == 'ELINSWVESQTNGIIR' or True:
            #    if modified_sequence not in peptides:
            #        peptides[modified_sequence] = {
            #            'id':peptide_count,
            #            'Protein IDs': "{}_{}".format(l['Proteins'], modified_sequence),
            #            'Reverse': '+' if 'REV_' in l['Proteins']  else '',
            #            'Potential contaminant': '+' if 'CON_' in l['Proteins'] else ''
            #        }
            #        peptide_count += 1
            #    l['Proteins'] = peptides[modified_sequence]['Protein IDs']
            #    l['Leading proteins'] = peptides[modified_sequence]['Protein IDs']
            #    l['Leading razor protein'] = peptides[modified_sequence]['Protein IDs']
            #    l['Protein group IDs'] = peptides[modified_sequence]['id']
            #if l['Type'] != 'MULTI-SECPEP':
            features.append(l)
    return features, peptides

def write_peptides(peptides, output_file):
    headers = ['id', 'Protein IDs', 'Reverse', 'Potential contaminant']
    with open(output_file, 'w') as w:
        r = DictWriter(w, fieldnames = headers, delimiter = '\t')
        r.writeheader()
        for peptide in peptides.values():
            r.writerow(peptide)

def write_evidence(features, output_file):
    headers = list(features[0].keys())
    with open(output_file, 'w') as w:
        r = DictWriter(w, fieldnames = headers, delimiter = '\t')
        r.writeheader()
        for feature in features:
            r.writerow(feature)

def read_all_peptides(all_peptides_txt,files_to_ids,file_mapping,scan_to_rt):
    feature = {}
    considered_feature = defaultdict(int)
    file_idx_count = defaultdict(int)
    feature_bucket = defaultdict(lambda: [[] for _ in range(2000)])
    with open(all_peptides_txt) as f:
        r = DictReader(f, delimiter = '\t')
        for i,l in enumerate(r):
            file_num = files_to_ids.index(file_mapping[l['Raw file']].split('/')[-1])
            file_idx = file_idx_count[file_num]
            file_idx_count[file_num] += 1
            considered_feature[(file_idx,file_num)] = 0
            l['Idx'] = i
            l['IdxInFile'] = file_idx
            l['FileNum'] = file_num
            feature[(file_idx,file_num)] = l
            l['Retention time start'] = scan_to_rt[(l['Raw file'],int(l['Min scan number']))]
            l['Retention time end'] = scan_to_rt[(l['Raw file'],int(l['Max scan number']))]
            l['Raw file'] = file_mapping[l['Raw file']].split('/')[-1]
            l['Variant count'] = 0
            l['Type'] = 'MULTI'
            l['Sequence'] = ' '
            l['Modified sequence'] = ' '
            l['Length'] = ' '
            l['Modifications'] = ' '
#             print(l)
#             raise Exception
#             print(int(float(l['m/z'])))
#             print(int(float(l['Retention time'])))
            feature_bucket[l['Raw file']][int(float(l['m/z']))].append(l)
    return feature, feature_bucket

def read_matched_features(matched_features_txt, features):
    aligned_features = []
    count = 0
    mapped = set()
    # single = 0
    with open(matched_features_txt) as f:
        r = DictReader(f, delimiter = '\t')
        for m,l in enumerate(r):
            aligned_features.append([])
            raw_files = l['Raw files'].split(";")
#             if len(l['Multiplet ids'].split(";")) == 1:
#                 single += 1
            for i,n in enumerate(l['Multiplet ids'].split(";")):
                mapped.add((int(n),int(raw_files[i])))
                features[(int(n),int(raw_files[i]))]['Mapping group'] = m
                aligned_features[-1].append(features[(int(n),int(raw_files[i]))])
    for key,feature in features.items():
        if key not in mapped:
            feature['Mapping group'] = -1
            aligned_features.append([feature])
    return aligned_features
