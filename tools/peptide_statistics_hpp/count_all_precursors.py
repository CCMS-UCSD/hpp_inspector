import sys
from python_ms_utilities import mztab
import csv
from pathlib import Path
from collections import defaultdict

output_combined = sys.argv[1]
input_metaworkflow_files = sys.argv[2:]

def get_decoys_from_mztab(mztab_file):
    filenames = {}
    all_decoy_outputs = ''
    output = defaultdict(lambda: defaultdict(lambda: defaultdict(set)))
    with open(mztab_file) as f:
        nextline = f.readline()
        while(nextline[0:3] == 'MTD'):
            if 'ms_run' in nextline:
                header_line = nextline.rstrip().split("\t")
                ms_filename = header_line[1].replace("-location","")
                ms_filepath = header_line[2].replace("file://","f.")
                filenames[ms_filename] = ms_filepath
            nextline = f.readline()
        while(nextline[0:3] == 'PRH' or nextline[0:3] == 'PRT' or nextline[0:3] == 'COM' or nextline == '\n'):
            nextline = f.readline()
        headers = nextline.rstrip().split('\t')
        mztab_dict = csv.DictReader(f, fieldnames = headers, delimiter = '\t')
        for row in list(mztab_dict):
            is_decoy = 'Target' if row['opt_global_cv_MS:1002217_decoy_peptide'] == 'null' else 'Decoy'
            try:
                is_canonical = 'Canonical' if float(row['opt_global_TopCanonicalProtQValue']) <= 0.01 else 'IsoformTrEMBL'
            except:
                is_canonical = 'IsoformTrEMBL'
            protein = row['accession']
            sequence = row['sequence']
            modifications = row['modifications']
            charge = row['charge']
            psm_id = row['PSM_ID']
            output['PSM'][is_canonical][is_decoy].add(psm_id)
            output['Precursor'][is_canonical][is_decoy].add((sequence,modifications,charge))
            output['Sequence'][is_canonical][is_decoy].add(sequence)
            output['Protein'][is_canonical][is_decoy].add(protein)
    counts = {}
    for fdr_kind, fdr_outputs in output.items():
        for is_canonical, fdr_outputs2 in fdr_outputs.items():
            for is_decoy, output_set in fdr_outputs2.items():
                counts["{}-{}-{}".format(is_decoy, is_canonical, fdr_kind)] = len(output_set)
    return counts

with open(output_combined,'w') as f:
    fieldnames = [
        'Job',
        'Description',
        'mzTab',
        'Decoy-Canonical-PSM',
        'Target-Canonical-PSM',
        'Decoy-IsoformTrEMBL-PSM',
        'Target-IsoformTrEMBL-PSM',
        'Decoy-Canonical-Precursor',
        'Target-Canonical-Precursor',
        'Decoy-IsoformTrEMBL-Precursor',
        'Target-IsoformTrEMBL-Precursor',
        'Decoy-Canonical-Sequence',
        'Target-Canonical-Sequence',
        'Decoy-IsoformTrEMBL-Sequence',
        'Target-IsoformTrEMBL-Sequence',
        'Decoy-Canonical-Protein',
        'Target-Canonical-Protein',
        'Decoy-IsoformTrEMBL-Protein',
        'Target-IsoformTrEMBL-Protein'
    ]
    w = csv.DictWriter(f, delimiter = '\t', fieldnames = fieldnames)
    for metaworkflow_file in input_metaworkflow_files:
        print(metaworkflow_file)
    for metaworkflow_file in input_metaworkflow_files:
        with open(metaworkflow_file) as metaworkflow_f:
            metaworkflow_lines = csv.DictReader(metaworkflow_f, delimiter='\t')
            for workflow_line in metaworkflow_lines:
                users = [workflow_line.get('User','batch'),'jswertz']
                for user in users:
                    try:
                        job = workflow_line['Amb_ID']
                        desc = workflow_line['Amb_desc']
                        mztab_path = '/data/ccms-data/tasks/{}/{}/mzTab/Top1_results.mzTab'.format(user,job)
                        output_dict = get_decoys_from_mztab(mztab_path)
                        print("{} ({})".format(desc,job))
                        print("Top1_results.mzTab")
                        for category,count in output_dict.items():
                            print('\t{}: {}'.format(category,count))
                        print('\n')
                        mztab_lib_path = '/data/ccms-data/tasks/{}/{}/mzTab_lib_construction/Top1_library_construction_results.mzTab'.format(user,job)
                        output_dict_lib = get_decoys_from_mztab(mztab_lib_path)
                        print("Top1_library_construction_results.mzTab")
                        for category,count in output_dict_lib.items():
                            print('\t{}: {}'.format(category,count))
                        print('\n')
                        output_dict['Job'] = job
                        output_dict_lib['Job'] = job
                        output_dict['Description'] = desc
                        output_dict_lib['Description'] = desc
                        output_dict['mzTab'] = 'Top1_results.mzTab'
                        output_dict_lib['mzTab'] = 'Top1_library_construction_results.mzTab'
                        w.writerow(output_dict)
                        w.writerow(output_dict_lib)
                    except:
                        pass
