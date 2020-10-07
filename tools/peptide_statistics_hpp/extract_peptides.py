import argparse
import sys
import csv
import glob
import requests
import urllib.parse
import os
from quantlib import mztab, proteosafe

def arguments():
    parser = argparse.ArgumentParser(description='mzTab to list of peptides')
    parser.add_argument('-r','--params', type = str, help='Input Parameters')
    parser.add_argument('-i','--mztab', type = str, help='Input Identifications')
    parser.add_argument('-p','--peptidelist', type = str, help='Output Peptides')
    parser.add_argument('-l','--output_library_folder', type = str, help='Output Library Folder')

    if len(sys.argv) < 4:
        parser.print_help()
        sys.exit(1)
    return parser.parse_args()

def convert_msv_to_task(msv):

    r = requests.get('https://massive.ucsd.edu/ProteoSAFe/FindDatasets?query={}'.format(msv))
    parsed = urllib.parse.urlparse(urllib.parse.unquote(r.url).split("url=")[-1])
    return urllib.parse.parse_qs(parsed.query)['task'][0]

def outputs_for_protein_explorer(psms, representatives_input_name, representatives_filtered_name, provenance_spectra_name):

    seen_representatives = set()

    representatives_input_fields = [
       'peptide', 'charge', 'length', 'mz', 'score', 'originalfile_filename',
       'originalfile_scan'
    ]

    representatives_filtered_fields = [
       'peptide', 'charge'
    ]

    provenance_spectra_fields = [
        'filename', 'annotation', 'charge', 'workflow', 'proteosafe_task',
        'proteosafe_moddecode_task', 'search_url',
        'scan', 'modifications'
    ]

    msv_to_task = {}

    with open(representatives_input_name, 'w') as representatives_input_file, open(representatives_filtered_name, 'w') as representatives_filtered_file, open(provenance_spectra_name, 'w') as provenance_spectra_file:

        representatives_input = csv.DictWriter(representatives_input_file, delimiter = '\t', fieldnames=representatives_input_fields)
        representatives_filtered = csv.DictWriter(representatives_filtered_file, delimiter = '\t', fieldnames=representatives_filtered_fields)
        provenance_spectra = csv.DictWriter(provenance_spectra_file, delimiter = '\t', fieldnames=provenance_spectra_fields)

        representatives_input.writeheader()
        representatives_filtered.writeheader()
        provenance_spectra.writeheader()

        for ((filename,scan),psm_list) in psms.items():
            psm = psm_list[0]
            peptide = psm.sequence
            charge = psm.charge
            msv = filename.split('/')[0].replace('f.','')
            if msv in msv_to_task:
                dataset_task = msv_to_task[msv]
            else:
                dataset_task = convert_msv_to_task(msv)
                msv_to_task[msv] = dataset_task
            provenance_spectra.writerow({
                'filename': filename[2:],
                'annotation': peptide,
                'charge': charge,
                'workflow': psm.kind,
                'proteosafe_task': dataset_task,
                'proteosafe_moddecode_task':'',
                'search_url': 'https://massive.ucsd.edu/ProteoSAFe/result.jsp?task={}&view=group_by_spectrum&file={}'.format(dataset_task, psm.search_filename),
                'scan':scan,
                'modifications':psm.modifications
            })
            if (peptide,charge) not in seen_representatives:
                seen_representatives.add((peptide,charge))
                representatives_input.writerow({
                    'peptide':peptide,
                    'charge':charge,
                    'length':len([p for p in peptide if p.isalpha()]),
                    'mz':psm.parent_mass,
                    'score':psm.score,
                    'originalfile_filename':filename[2:],
                    'originalfile_scan':scan
                })
                representatives_filtered.writerow({
                    'peptide':peptide,
                    'charge':charge
                })



def main():
    args = arguments()
    ids = {}
    file_mapping = proteosafe.read_params(args.params, full_paths = True)
    print(file_mapping)
    if args.mztab != "" and args.mztab != None:
        for mztab_file in glob.glob(args.mztab + '/*'):
            try:
                ids = mztab.read(mztab_file, ids, file_mapping[mztab_file.split('/')[-1]])
            except:
                ids = mztab.read_lib(mztab_file, ids)

    peptides = set([
        ''.join([aa for aa in id[0].sequence if aa.isalpha()])
        for filescan,id in
        ids.items()
        if not ('j10022_RNF14' in filescan[0] or 'j12271_TMPRSS3' in filescan[0] or 'j13961_ADAM32' in filescan[0] or 'j9837_ADCYAP1' in filescan[0] or 'q8920_ACCN4' in filescan[0])
    ])

    with open(args.peptidelist,'w') as w:
        w.write('Peptide\n')
        for peptide in peptides:
            w.write('{}\n'.format(peptide))

    representatives_input_name = os.path.join(args.output_library_folder,"representatives_input.tsv")
    representatives_filtered_name = os.path.join(args.output_library_folder,"representatives_filtered.tsv")
    provenance_spectra_name = os.path.join(args.output_library_folder,"provenance_spectra.tsv")

    outputs_for_protein_explorer(ids, representatives_input_name, representatives_filtered_name, provenance_spectra_name)


if __name__ == '__main__':
    main()
