import argparse
import sys
import csv
from pathlib import Path
from collections import defaultdict
import read_mappings
from python_ms_utilities import mapping, resources

def arguments():
    parser = argparse.ArgumentParser(description='mzTab to list of peptides')
    parser.add_argument('-p','--params', type = str, help='Input Parameters')
    parser.add_argument('-c','--comparisons', type = Path, help='Comparison Jobs')
    parser.add_argument('-f','--proteome_fasta', type = Path, help='FASTA File')
    parser.add_argument('-b','--backup_kb_pep', type = Path, help='Backup KB Peptides')
    parser.add_argument('-t','--use_job_level_thresholds', type = bool, help='Use job level thresholds',default=True)
    parser.add_argument('-k','--kb_pep', type = Path, help='Output KB Peptides')
    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit(1)
    return parser.parse_args()

def read_coverage_folder(input_folder,proteome, job_level_threshold):
    
    protein_mapping_out_combined = {}

    protein_mapping_out_per_file = []
    protein_hpp_fdr_per_file = []
    protein_hint_fdr_per_file = []

    all_proteins = []

    for protein_coverage_file in input_folder.glob('*'):
        protein_mapping_out,_,_,_,protein_hpp_fdr,protein_hint_fdr = read_mappings.read_protein_coverage(protein_coverage_file,set(),proteome,True,False,job_level_threshold)
        protein_mapping_out_per_file.append(protein_mapping_out)
        protein_hpp_fdr_per_file.append(protein_hpp_fdr)
        protein_hint_fdr_per_file.append(protein_hint_fdr)
        all_proteins.extend(list(protein_mapping_out.keys()))

    all_proteins = set(all_proteins)

    for protein in all_proteins:
        file_to_consider_per_protein = -1
        best_hpp_fdr_per_protein = 1
        best_hint_fdr_per_protein = 1
        files_pass_hpp_fdr = [(i,fdr.get(protein,1)) for i,fdr in enumerate(protein_hpp_fdr_per_file) if float(fdr.get(protein,1)) <= 0.01]
        files_pass_hint_fdr = [(i,fdr.get(protein,1)) for i,fdr in enumerate(protein_hint_fdr_per_file) if float(fdr.get(protein,1)) <= 0.01]
        if len(files_pass_hpp_fdr) > 0:
            best_hpp_fdr_per_protein = min([f[1] for f in files_pass_hpp_fdr])
            file_to_consider_per_protein = [f[0] for f in files_pass_hpp_fdr if f[1] == best_hpp_fdr_per_protein][0]
        elif len(files_pass_hint_fdr) > 0:
            #if only hints, just keep lowest FDR set
            best_hint_fdr_per_protein = min([f[1] for f in files_pass_hint_fdr])
            file_to_consider_per_protein = [f[0] for f in files_pass_hint_fdr if f[1] == best_hint_fdr_per_protein][0]
        if file_to_consider_per_protein >= 0:
            protein_mapping_out_combined[protein] = protein_mapping_out_per_file[file_to_consider_per_protein][protein]

    return protein_mapping_out_combined

def main():
    args = arguments()

    header = ['protein','aa_start','aa_end','demodified','synthetic_cosine','all_protein_fdr','hpp_protein_fdr','is_hpp']

    if args.comparisons:

        try:

            proteome = mapping.add_decoys(mapping.read_uniprot(args.proteome_fasta))

            with open(args.kb_pep, 'w') as w:
                r = csv.DictWriter(w, delimiter = '\t', fieldnames = header)
                r.writeheader()
                protein_mapping_out = read_coverage_folder(args.comparisons, proteome, args.use_job_level_thresholds)
                for protein, peptide_mappings in protein_mapping_out.items():
                    for peptide, mappings in peptide_mappings.items():
                        for (start, end, cosine, all_protein_fdr, hpp_protein_fdr, is_hpp) in mappings:
                            r.writerow({
                                'protein':protein,
                                'aa_start':start,
                                'aa_end':end,
                                'demodified':peptide,
                                'synthetic_cosine':cosine,
                                'all_protein_fdr':all_protein_fdr,
                                'hpp_protein_fdr':hpp_protein_fdr,
                                'is_hpp':is_hpp
                            })
        except:
            with open(args.kb_pep, 'w') as w:
                r = csv.DictWriter(w, delimiter = '\t', fieldnames = header)
                r.writeheader()

                in_vivo = set()
                synthetic = set()

                for kb_input in args.comparisons.glob('*'):
                    with open(kb_input) as f:
                        kb_rs = csv.DictReader(f, delimiter = '\t')
                        for kb_row in kb_rs:
                            if kb_row['library'] == '2':
                                in_vivo.add((
                                    kb_row['protein'],
                                    kb_row['aa_start'],
                                    kb_row['aa_end'],
                                    kb_row['demodified']
                                ))
                            else:
                                synthetic.add((
                                    kb_row['protein'],
                                    kb_row['aa_start'],
                                    kb_row['aa_end'],
                                    kb_row['demodified']
                                ))
                for in_vivo_seq in in_vivo:
                    r.writerow({
                        'protein':in_vivo_seq[0],
                        'aa_start':in_vivo_seq[1],
                        'aa_end':in_vivo_seq[2],
                        'demodified':in_vivo_seq[3],
                        'synthetic_cosine': '0' if in_vivo_seq in synthetic else '-1',
                        'all_protein_fdr': '1',
                        'hpp_protein_fdr': '0',
                        'is_hpp':'True'
                    })
    else:
        args.kb_pep.write_text(args.backup_kb_pep.read_text())

    # try:
    #
    #     PAGE_SIZE = 500000;
    #     OFFSET = 0;
    #
    #     peptides_to_add = requests.get('http://ccms-internal.ucsd.edu/ProteoSAFe/PeptideStatisticsServlet?pageSize={}&offset={}'.format(PAGE_SIZE,OFFSET)).json()
    #     peptides = peptides_to_add
    #     while (len(peptides_to_add) == PAGE_SIZE):
    #         OFFSET += PAGE_SIZE;
    #         peptides_to_add = requests.get('http://ccms-internal.ucsd.edu/ProteoSAFe/PeptideStatisticsServlet?pageSize={}&offset={}'.format(PAGE_SIZE,OFFSET)).json()
    #         peptides.extend(peptides_to_add)
    #
    #     header = ['library','protein','aa_start','aa_end','demodified','charge', 'sequence']
    #     with open(args.kb_pep, 'w') as w:
    #         r = csv.DictWriter(w, delimiter = '\t', fieldnames = header)
    #         r.writeheader()
    #         for peptide in peptides:
    #             r.writerow(peptide)
    #
    # except:



if __name__ == '__main__':
    main()
