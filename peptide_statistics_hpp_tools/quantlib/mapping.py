from collections import defaultdict
from quantlib import mass
from csv import DictWriter

def convert_features(features):
    features_converted = {}
    for file, feature_bins in features.items():
        for feature_bin in feature_bins:
            for feature in feature_bin:
                features_converted[(feature['IdxInFile'], feature['FileNum'])] = feature
    return features_converted

def write_features(features, feature_output_file, include_noid = False):
    feature_header = ['id','Raw file', 'Type', 'Mapping group', 'Charge', 'm/z', 'Mass', 'Uncalibrated m/z', 'Resolution', 'Number of data points', 'Number of scans', 'Number of isotopic peaks', 'PIF', 'Mass fractional part', 'Mass deficit', 'Mass precision [ppm]', 'Max intensity m/z 0', 'Retention time', 'Retention length', 'Retention length (FWHM)', 'Retention time start', 'Retention time end', 'Min scan number', 'Max scan number', 'Identified', 'MS/MS IDs', 'Variant count', 'Sequence', 'Length', 'Modifications', 'Modified sequence', 'Proteins', 'Protein group IDs','Reverse','Potential contaminant', 'Score', 'Intensity', 'Intensities', 'Isotope pattern', 'MS/MS Count', 'MSMS Scan Numbers', 'MSMS Isotope Indices', 'Idx', 'IdxInFile', 'FileNum', 'MSMS Mapped Isotope Indices', 'MSMS Mapped Scan Numbers', 'MSMS Charge', 'Map to Different File']
    with open(feature_output_file, "w") as w:
        r = DictWriter(w, feature_header, delimiter = '\t')
        r.writeheader()
        for i,feature in enumerate(features):
            if include_noid or feature['Sequence'] != ' ':
                feature['id'] = i
                r.writerow(feature)

def msms_to_features(features, ids, mz_tol, rt_tol, use_all_mappings):
    msms_to_rtmz = {}
    features_w_msms = 0
    feature_to_msms = defaultdict(list)
    isotopes = []
    errors= []
    successes= []
    comparisons = 0
    count = 0
    proteins = {}
    protein_count = 0
    ms2_closest_to_feature = defaultdict(lambda: 500)
    for p,(spec,peptides) in enumerate(ids.items()):
        best_scans = 100
        best_mz = 0
        best_feature = None
        best_iso = 100
        best_run = None
        error = None
        success = None
        if (p+1) % 1000 == 0:
            print("{}/{}".format(p+1,len(ids)))
        for peptide in peptides:
            for run in features:
                ms2_from_run = run == spec[0]
                if ms2_from_run or use_all_mappings:
                    pep_mz = mass.theoretical_mz(peptide.sequence,peptide.charge)
                    feature_buckets = []
                    for i,fs in enumerate(features[run][int(pep_mz) - 2:int(pep_mz) + 3]):
                        for k,feature in enumerate(fs):
                            comparisons += 1
                            feature_mz = float(feature['m/z'])
                            feature_rt = float(feature['Retention time'])
                            feature_charge = int(feature['Charge'])
                            if feature_charge == peptide.charge:
                                rt = peptide.rt
                                feature_min_rt = float(feature['Retention time start'])*60
                                feature_max_rt = float(feature['Retention time end'])*60
                                if rt:
                                    if rt <= feature_max_rt and rt >= feature_min_rt:
                                        percent_off = 0
                                    else:
                                        percent_off = min(abs(rt-feature_max_rt),abs(rt-feature_min_rt))/(feature_max_rt-feature_min_rt)
                                else:
                                    feature_min_scan = int(feature['Min scan number'])
                                    feature_max_scan = int(feature['Max scan number'])
                                    if spec[1] <= feature_max_scan and spec[1] >= feature_min_scan:
                                        percent_off = 0
                                    else:
                                        percent_off = min(abs(spec[1]-feature_min_scan),abs(spec[1]-feature_max_scan))/(feature_max_scan-feature_min_scan)
                                ppm_difference = min(list(enumerate([
        #                             abs(feature_mz - pep_mz - 4/feature_charge),
        #                             abs(feature_mz - pep_mz - 3/feature_charge),
        #                             abs(feature_mz - pep_mz - 2/feature_charge),
        #                             abs(feature_mz - pep_mz - 1/feature_charge),
                                    abs(feature_mz - pep_mz),
                                    abs(feature_mz - pep_mz + 1/feature_charge),
                                    abs(feature_mz - pep_mz + 2/feature_charge),
                                    abs(feature_mz - pep_mz + 3/feature_charge),
                                    abs(feature_mz - pep_mz + 4/feature_charge)
                                ])), key = lambda x: x[1])
                                isotope = ppm_difference[0]
                                ppm_difference = (ppm_difference[1]/pep_mz)*1e6
        #                         if features[spec[0]][i + int(pep_mz) - 2][k]['Idx'] not in feature_to_msms:
                                if ppm_difference < mz_tol and percent_off <= rt_tol and isotope == 0:
                                        best_scans = percent_off
                                        best_mz = feature_mz
                                        best_feature = (i + int(pep_mz) - 2,k)
                                        best_iso = 0
                                        success = (percent_off, ppm_difference, isotope, spec, peptide, feature, ms2_from_run)
                                        best_run = run
        #                         elif ppm_difference < mz_tol and percent_off <= rt_tol:
        #                             if percent_off < best_scans:
        #                                 best_scans = percent_off
        #                                 best_mz = feature_mz
        #                                 best_feature = (i + int(pep_mz) - 2,k)
        #                                 best_iso = isotope
        #                                 success = (percent_off, ppm_difference, isotope, spec, peptide, feature, peptide[3])
                                if ppm_difference < 1000 and percent_off < .5:
                                    if error:
                                        closest_ppm = min(ppm_difference,error[1])
                                        if closest_ppm == ppm_difference:
                                            error = (percent_off, ppm_difference, isotope, spec, peptide, feature, ms2_from_run)
                                    else:
                                        error = (percent_off, ppm_difference, isotope, spec, peptide, feature, ms2_from_run)
        if best_feature:
            count += 1
#             features[spec[0]][best_feature[0]][best_feature[1]]['Peptide'] = peptide.sequence
#             features[spec[0]][best_feature[0]][best_feature[1]]['Charge'] = peptide.charge
            msms_to_rtmz[spec] = features[best_run][best_feature[0]][best_feature[1]]
            feature_to_msms[features[best_run][best_feature[0]][best_feature[1]]['Idx']].append((spec, best_iso))
            modified_sequence = success[4][0]
            sequence = "".join([s for s in modified_sequence if s.isalpha()])
            # feature_assigned_proteins = success[4].protein + "_" + modified_sequence
            feature_assigned_proteins = success[4].protein
            if success[1] < ms2_closest_to_feature[features[best_run][best_feature[0]][best_feature[1]]['Idx']]:
                if feature_assigned_proteins not in proteins:
                    proteins[feature_assigned_proteins] = {
                        'id':protein_count,
                        'Protein IDs': feature_assigned_proteins,
                        'Reverse': '+' if 'REV_' in feature_assigned_proteins  else '',
                        'Potential contaminant': '+' if 'CON_' in feature_assigned_proteins else ''
                    }
                    protein_count += 1
                ms2_closest_to_feature[features[best_run][best_feature[0]][best_feature[1]]['Idx']] = success[1]
                features[best_run][best_feature[0]][best_feature[1]]['Protein group IDs'] = proteins[feature_assigned_proteins]['id']
                features[best_run][best_feature[0]][best_feature[1]]['Proteins'] = feature_assigned_proteins
                features[best_run][best_feature[0]][best_feature[1]]['Reverse'] = proteins[feature_assigned_proteins]['Reverse']
                features[best_run][best_feature[0]][best_feature[1]]['Potential contaminant'] = proteins[feature_assigned_proteins]['Potential contaminant']
            modifications = success[4].modifications
            kind = success[4].kind
            charge = success[4].charge
            msms_scan = spec[1]
            isotope_indicies = success[2]
            id_conflict = False
            if features[best_run][best_feature[0]][best_feature[1]]['Sequence'] != " " and features[best_run][best_feature[0]][best_feature[1]]['Sequence'] != sequence:
                id_conflict = True

            # features[spec[0]][best_feature[0]][best_feature[1]]['Proteins'] = proteins
            if not id_conflict:
                features[best_run][best_feature[0]][best_feature[1]]['Variant count'] = 1
                features[best_run][best_feature[0]][best_feature[1]]['Type'] = kind
                features[best_run][best_feature[0]][best_feature[1]]['Sequence'] = sequence
                features[best_run][best_feature[0]][best_feature[1]]['Modified sequence'] = modified_sequence
                features[best_run][best_feature[0]][best_feature[1]]['Length'] = str(len(sequence))
                features[best_run][best_feature[0]][best_feature[1]]['Modifications'] = modifications
            else:
                features[best_run][best_feature[0]][best_feature[1]]['Variant count'] += 1
                features[best_run][best_feature[0]][best_feature[1]]['Type'] = "IDCONFLICT_" + features[best_run][best_feature[0]][best_feature[1]]['Type'].replace("MULTI_","").replace("IDCONFLICT_","")
                features[best_run][best_feature[0]][best_feature[1]]['Sequence'] += ";" + sequence
                features[best_run][best_feature[0]][best_feature[1]]['Modified sequence'] += ";" + modified_sequence
                features[best_run][best_feature[0]][best_feature[1]]['Length'] += ";" + str(len(sequence))
                features[best_run][best_feature[0]][best_feature[1]]['Modifications'] += ";" + modifications

            if 'MSMS Mapped Scan Numbers' not in features[best_run][best_feature[0]][best_feature[1]]:
                features[best_run][best_feature[0]][best_feature[1]]['MSMS Mapped Scan Numbers'] = str(msms_scan)
                features[best_run][best_feature[0]][best_feature[1]]['MSMS Mapped Isotope Indices'] = str(isotope_indicies)
                features[best_run][best_feature[0]][best_feature[1]]['MSMS Charge'] = str(charge)
                features[best_run][best_feature[0]][best_feature[1]]['Map to Different File'] = ms2_from_run
            elif features[best_run][best_feature[0]][best_feature[1]]['MSMS Mapped Scan Numbers'] != "":
                features[best_run][best_feature[0]][best_feature[1]]['MSMS Mapped Scan Numbers'] += ";" + str(msms_scan)
                features[best_run][best_feature[0]][best_feature[1]]['MSMS Mapped Isotope Indices'] += ";" + str(isotope_indicies)
                features[best_run][best_feature[0]][best_feature[1]]['MSMS Charge'] += ";" + str(charge)
                features[best_run][best_feature[0]][best_feature[1]]['Map to Different File'] = features[best_run][best_feature[0]][best_feature[1]]['Map to Different File'] or ms2_from_run
            else:
                features[best_run][best_feature[0]][best_feature[1]]['MSMS Mapped Scan Numbers'] = str(msms_scan)
                features[best_run][best_feature[0]][best_feature[1]]['MSMS Mapped Isotope Indices'] = str(isotope_indicies)
                features[best_run][best_feature[0]][best_feature[1]]['MSMS Charge'] = str(charge)
                features[best_run][best_feature[0]][best_feature[1]]['Map to Different File'] = ms2_from_run
            features[best_run][best_feature[0]][best_feature[1]]['MS/MS Count'] = str(int(features[best_run][best_feature[0]][best_feature[1]]['MS/MS Count']) + 1)
            successes.append(success)
        elif error:
            errors.append(error)
    print(count)
    return features, proteins

def propogate_over_matched_features(matched_features):
    output_features = []
    for features in matched_features:
        id_feature = None
        is_misaligned = False
        for feature in features:
            feature_sequence = feature['Modified sequence'].replace('I','L')
            if id_feature and feature_sequence != ' ' and id_feature['Modified sequence'].replace('I','L') != feature_sequence:
                is_misaligned = True
            if feature_sequence != ' ' and  feature_sequence != '':
                id_feature = feature
        for feature in features:
            if feature['Sequence'] == ' ' and id_feature:
                feature['Sequence'] = id_feature['Sequence']
                feature['Modified sequence'] = id_feature['Modified sequence']
                feature['Length'] = id_feature['Length']
                feature['Modifications'] = id_feature['Modifications']
                feature['Proteins'] = id_feature['Proteins']
                feature['Protein group IDs'] = id_feature['Protein group IDs']
                feature['Reverse'] = id_feature['Reverse']
                feature['Potential contaminant'] = id_feature['Potential contaminant']
                if is_misaligned:
                    feature['Type'] = 'MAPCONFLICT_' + id_feature['Type']
                else:
                    feature['Type'] = 'MAPPED_' + id_feature['Type']
            output_features.append(feature)
    return output_features
