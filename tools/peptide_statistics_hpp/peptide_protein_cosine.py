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
from pyteomics import mzml
from pathlib import Path
import pickle

def arguments():
    parser = argparse.ArgumentParser(description='mzTab to list of peptides')
    parser.add_argument('-k','--kb_pep', type = str, help='Peptides from KB')
    parser.add_argument('-n','--input_proteins', type = Path, help='Input Proteins')
    parser.add_argument('-x','--input_psms', type = Path, help='Input PSMs')
    parser.add_argument('-e','--output_psms', type = Path, help='Output PSMs')
    parser.add_argument('-p','--output_peptides', type = Path, help='Output Peptides')
    parser.add_argument('-r','--output_proteins', type = Path, help='Output Proteins')
    parser.add_argument('-c','--protein_coverage', type = str, help='Added Protein Coverage')
    parser.add_argument('-t','--cosine_cutoff', type = float, help='Cosine Cutoff')
    parser.add_argument('-l','--explained_intensity_cutoff', type = float, help='Explained Intensity Cutoff')
    parser.add_argument('-z','--nextprot_pe', type = str, help='NextProt PEs')

    if len(sys.argv) < 4:
        parser.print_help()
        sys.exit(1)
    return parser.parse_args()

def find_overlap(existing_peptides, new_peptides, protein_length, protein_pe, name):
    kb_pos = set([sorted(positions, key = lambda x: x[0])[0] for positions in existing_peptides.values()])
    added_pos = set([sorted(positions, key = lambda x: x[0])[0] for positions in new_peptides.values()])

    kb_peptides = set(["{} ({}-{})".format(peptide,sorted(positions, key = lambda x: x[0])[0][0],sorted(positions, key = lambda x: x[0])[0][1]) for peptide, positions in existing_peptides.items()])
    added_peptides = set(["{} ({}-{})".format(peptide,sorted(positions, key = lambda x: x[0])[0][0],sorted(positions, key = lambda x: x[0])[0][1]) for peptide, positions in new_peptides.items()])

    novel_peptides = added_peptides.difference(kb_peptides)
    supporting_peptides = added_peptides.intersection(kb_peptides)

    nonoverlapping_kb_peptides = hupo_nonoverlapping(kb_pos.difference(added_pos))
    nonoverlapping_added_peptides = hupo_nonoverlapping(added_pos.difference(kb_pos))
    nonoverlapping_intersection = hupo_nonoverlapping(kb_pos.union(added_pos))
    nonoverlapping_all_kb = hupo_nonoverlapping(kb_pos)

    return {
        'supporting_peptides'+name:','.join(supporting_peptides) if len(supporting_peptides) >= 1 else ' ',
        'novel_peptides'+name:','.join(novel_peptides) if len(novel_peptides) >= 1 else ' ',
        'kb_hpp'+name:nonoverlapping_all_kb,
        'new_hpp'+name:nonoverlapping_added_peptides,
        'combined_hpp'+name:nonoverlapping_intersection,
        'added_kb_coverage'+name:find_coverage(list(existing_peptides.values()) + list(new_peptides.values()),protein_length),
        'promoted'+name:'Yes' if (nonoverlapping_all_kb < 2 and nonoverlapping_intersection >= 2) else 'No'
    }, [s.split(' ')[0] for s in supporting_peptides]

def find_coverage(pos, length):
    coverage = [False for _ in range(length)]
    for p in pos:
        for (start,end) in p:
            for i in range(start,end+1):
                if i <= len(coverage):
                    # print(protein)
                    # print(i)
                    # print(len(coverage))
                    coverage[i-1] = True
    return len([c for c in coverage if c])

def hupo_nonoverlapping(segments, just_noncontained = True):
    if len(segments) == 0:
        return 0
    max_peps = [0 for _ in range(len(segments))]
    max_score_at_pep = [0 for _ in range(len(segments))]
    segments = [(int(s[0]),int(s[1])) for s in segments]
    ordered_segments = sorted(segments, key = lambda x: (int(x[1]),int(x[0])))

    for i in range(len(segments)):
        max_peps_j = 0
        max_score_at_pep_j = 0

        j = None
        for k in range(i-1,-1,-1):
            if (ordered_segments[k][0] <= ordered_segments[i][0] and
                ordered_segments[k][1] < ordered_segments[i][1] and
                ordered_segments[k][1] >= ordered_segments[i][0]):
                if j:
                    if (max_peps[k] > max_peps[j]):
                        j = k
                else:
                    j = k

        if j:
            max_peps_j = max_peps[j]

        j = None
        for k in range(i-1,-1,-1):
            if (ordered_segments[k][1] <= ordered_segments[i][0]):
                j = k
                break

        if j is not None:
            max_score_at_pep_j = max_score_at_pep[j]

        max_peps[i] = 1 + max(max_peps_j,max_score_at_pep_j)

        max_score_at_pep_j = 0

        j = None
        for k in range(i-1,-1,-1):
            if (ordered_segments[k][1] < ordered_segments[i][1] and
                ordered_segments[k][1] >= ordered_segments[i][0]):
                if j:
                    if (max_score_at_pep[k] > max_score_at_pep[j]):
                        j = k
                else:
                    j = k

        if j is not None:
            max_score_at_pep_j = max_score_at_pep[j]

        max_score_at_pep[i] = max(max_score_at_pep_j,max_peps[i])

    return max_score_at_pep[-1]

def main():
    args = arguments()

    representative_per_precursor = {}
    protein_length = {}

    added_proteins = defaultdict(lambda: defaultdict(list))
    added_proteins_matching_synthetic = defaultdict(lambda: defaultdict(list))
    added_proteins_explained_intensity = defaultdict(lambda: defaultdict(list))

    # note that this has a limit of 100000000 so with multi-species in ProteinExplorer this will likely need to change
    url = "http://massive.ucsd.edu/ProteoSAFe/ProteinLibraryServlet?task=protein_explorer_proteins&file=&pageSize=100000000&offset=0&query=%2523%257B%2522unroll%2522%253A%2522no%2522%252C%2522include_synthetics%2522%253A%2522yes%2522%252C%2522datasets%2522%253A%2522%2522%252C%2522accession_input%2522%253A%2522%2522%257D&query_type=representative&_=1560375217897"
    proteins = requests.get(url).json()['row_data']
    for protein in proteins:
        # nextprot_pe[protein['accession']] = int(protein['proteinexistance'])
        protein_length[protein['accession']] = int(protein['length'])
    kb_proteins = defaultdict(lambda: defaultdict(list))
    added_proteins = defaultdict(lambda: defaultdict(list))
    all_proteins = defaultdict(lambda: defaultdict(list))
    kb_proteins_w_synthetic = defaultdict(lambda: defaultdict(list))
    added_proteins_w_synthetic = defaultdict(lambda: defaultdict(list))
    has_synthetic = set()

    with open(args.kb_pep) as f:
        r = csv.DictReader(f, delimiter='\t')
        for l in r:
            if int(l['library']) == 2:
                kb_proteins[l['protein']][l['demodified'].replace('I','L')].append((int(l['aa_start']),int(l['aa_end'])))
            if int(l['library']) == 3 or int(l['library']) == 4:
                has_synthetic.add(l['demodified'].replace('I','L'))

    for protein in kb_proteins:
        for seq in kb_proteins[protein].keys():
            if seq in has_synthetic:
                kb_proteins_w_synthetic[protein][seq] = kb_proteins[protein][seq]

    ms_existence = {}
    nextprot_pe = {}
    peptide_to_protein = {}
    peptide_to_protein_added = {}


    with open(args.nextprot_pe) as f:
        r = csv.DictReader(f, delimiter='\t')
        for l in r:
            ms_existence[l['protein'].replace('NX_','')] = l['ms_evidence']
            nextprot_pe[l['protein'].replace('NX_','')] = int(l['pe'])

    with open(args.protein_coverage) as f:
        r = csv.DictReader(f, delimiter='\t')
        for l in r:
            proteins_w_coords = l['proteins_w_coords']
            peptide = l['il_peptide']
            gene_unique = l['gene_unique'] == 'True'
            proteins = [
                [p.split(' ')[0]] + p.split(' ')[1].replace('(','').replace(')','').split('-')
                for p in
                proteins_w_coords.split(';')
            ]
            for protein in proteins:
                if len(protein) == 3:
                    accession = protein[0].split('|')[1]
                    aa_start = int(protein[1])
                    aa_end = int(protein[2])
                    all_proteins[accession][peptide].append((aa_start,aa_end))
                    if accession in nextprot_pe and (aa_end - aa_start >= 9) and gene_unique:
                        added_proteins[accession][peptide].append((aa_start,aa_end))

    for protein in all_proteins:
        for peptide in all_proteins[protein].keys():
            peptide_to_protein[peptide.replace('I','L')] = protein

    for protein in added_proteins:
        for peptide in added_proteins[protein].keys():
            peptide_to_protein_added[peptide.replace('I','L')] = protein
            if peptide.replace('I','L') in has_synthetic:
                added_proteins_w_synthetic[protein][peptide] = added_proteins[protein][peptide]

    print(peptide_to_protein)

    with open(args.output_psms,'w') as w:
        header = ['protein','pe','ms_evidence','filename','scan','sequence','charge','usi','score','pass','type','parent_mass','synthetic_filename','synthetic_scan','synthetic_usi','cosine','synthetic_match','explained_intensity','hpp_match']

        o = csv.DictWriter(w, delimiter='\t',fieldnames = header)
        o.writeheader()
        for input_psm in args.input_psms.glob('*'):
            with open(input_psm) as f:
                r = csv.DictReader(f, delimiter='\t')
                for l in r:
                    il_peptide = ''.join([p.replace('I','L') for p in l['sequence'] if p.isalpha()])
                    protein = peptide_to_protein[il_peptide]
                    l.update({
                        'protein': protein,
                        'pe': nextprot_pe[protein],
                        'ms_evidence':ms_existence.get(protein,'no')
                    })
                    if il_peptide in peptide_to_protein_added:
                        l['hpp_match'] = 'Yes'
                    else:
                        l['hpp_match'] = 'No'
                    o.writerow(l)
                    sequence, charge = l['sequence'],l['charge']
                    if (sequence, charge) in representative_per_precursor:
                        precursor_representative = representative_per_precursor[(sequence, charge)]
                        best_cosine = float(l['cosine']) > precursor_representative['cosine']
                        best_score = float(l['score']) > precursor_representative['score']
                        if best_cosine:
                            precursor_representative['cosine_filename'] = l['filename']
                            precursor_representative['cosine_scan'] = l['scan']
                            precursor_representative['cosine_usi'] = l['usi']
                            precursor_representative['synthetic_filename'] = l['synthetic_filename']
                            precursor_representative['synthetic_scan'] = l['synthetic_scan']
                            precursor_representative['synthetic_usi'] = l['synthetic_usi']
                            precursor_representative['cosine'] = float(l['cosine'])
                            precursor_representative['explained_intensity'] = float(l['explained_intensity'])
                        if best_score:
                            precursor_representative['database_filename'] = l['filename']
                            precursor_representative['database_scan'] = l['scan']
                            precursor_representative['database_usi'] = l['usi']
                            precursor_representative['score'] = float(l['score'])
                            precursor_representative['explained_intensity'] = float(l['explained_intensity'])
                        if (best_score and best_cosine):
                            precursor_representative['cosine_score_match'] = 'Yes'
                        else:
                            precursor_representative['cosine_score_match'] = 'No'

                    else:
                        l['database_filename'] = l['filename']
                        l['database_scan'] = l['scan']
                        l['database_usi'] = l['usi']
                        l['cosine_filename'] = l['filename']
                        l['cosine_scan'] = l['scan']
                        l['cosine_usi'] = l['usi']
                        l['explained_intensity'] = float(l['explained_intensity'])
                        l['cosine'] = float(l['cosine'])
                        l['score'] = float(l['score'])
                        l.pop('filename')
                        l.pop('scan')
                        l.pop('usi')
                        representative_per_precursor[(sequence, charge)] = l

    with open(args.output_peptides,'w') as w:
        header = ['protein','pe','ms_evidence','aa_total','database_filename','database_scan','database_usi','sequence','charge','score','pass','type','parent_mass','cosine_filename','cosine_scan','cosine_usi','synthetic_filename','synthetic_scan','synthetic_usi','cosine','synthetic_match','cosine_score_match','explained_intensity','hpp_match']
        r = csv.DictWriter(w, delimiter = '\t', fieldnames = header)
        r.writeheader()
        for (sequence, charge), best_psm in representative_per_precursor.items():
            sequence_il = ''.join([p.replace('I','L') for p in sequence if p.isalpha()])
            protein = peptide_to_protein[sequence_il]
            best_psm.update({
                'protein': protein,
                'pe': nextprot_pe[protein],
                'ms_evidence':ms_existence.get(protein,'no'),
                'aa_total':protein_length.get(protein,0)
            })
            if il_peptide in peptide_to_protein_added:
                l['hpp_match'] = 'Yes'
            else:
                l['hpp_match'] = 'No'
            r.writerow(best_psm)
            if sequence_il in added_proteins[best_psm['protein']]:
                if best_psm['explained_intensity'] > args.explained_intensity_cutoff:
                    if best_psm['cosine'] > args.cosine_cutoff:
                        added_proteins_matching_synthetic[best_psm['protein']][sequence_il] = added_proteins[best_psm['protein']][sequence_il]
                    added_proteins_explained_intensity[best_psm['protein']][sequence_il] = added_proteins[best_psm['protein']][sequence_il]

    with open(args.output_proteins, 'w') as fo:

        fieldnames = [
            'protein',
            'gene',
            'pe',
            'aa_total',
            'ms_evidence',
            'fdr',
            'supporting_peptides',
            'novel_peptides',
            'kb_hpp',
            'new_hpp',
            'combined_hpp',
            'added_kb_coverage',
            'supporting_peptides_w_synthetic',
            'novel_peptides_w_synthetic',
            'kb_hpp_w_synthetic',
            'new_hpp_w_synthetic',
            'combined_hpp_w_synthetic',
            'added_kb_coverage_w_synthetic',
            'supporting_peptides_w_synthetic_cosine',
            'novel_peptides_w_synthetic_cosine',
            'kb_hpp_w_synthetic_cosine',
            'new_hpp_w_synthetic_cosine',
            'combined_hpp_w_synthetic_cosine',
            'added_kb_coverage_w_synthetic_cosine',
            'promoted',
            'promoted_w_synthetic',
            'promoted_w_synthetic_cosine'
        ]

        w = csv.DictWriter(fo, delimiter = '\t', fieldnames = fieldnames)
        w.writeheader()

        for protein in added_proteins.keys():
            if (protein in nextprot_pe):
                protein_dict = {
                    'protein': protein,
                    'pe': nextprot_pe[protein],
                    'ms_evidence':ms_existence.get(protein,'no'),
                    'aa_total':protein_length.get(protein,0)
                }
                protein_dict.update(find_overlap(kb_proteins[protein],added_proteins_explained_intensity[protein_dict['protein']],int(protein_dict['aa_total']),int(protein_dict['pe']),'')[0])
                protein_dict.update(find_overlap(kb_proteins_w_synthetic[protein],added_proteins_matching_synthetic[protein_dict['protein']],int(protein_dict['aa_total']),int(protein_dict['pe']),'_w_synthetic')[0])
                protein_dict.update(find_overlap(kb_proteins_w_synthetic[protein],added_proteins_matching_synthetic[protein_dict['protein']],int(protein_dict['aa_total']),int(protein_dict['pe']),'_w_synthetic_cosine')[0])
                w.writerow(protein_dict)

if __name__ == '__main__':
    main()
