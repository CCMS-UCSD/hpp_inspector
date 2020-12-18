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
from itertools import chain

def arguments():
    parser = argparse.ArgumentParser(description='mzTab to list of peptides')
    parser.add_argument('-k','--kb_pep', type = Path, help='Peptides from KB')
    parser.add_argument('-s','--skip_kb', type = Path, help='Skip Peptides from KB')
    parser.add_argument('-x','--input_psms', type = Path, help='Input PSMs')
    parser.add_argument('-y','--input_psms_external', type = Path, help='Input PSMs (External)')
    parser.add_argument('-e','--output_psms', type = Path, help='Output PSMs')
    parser.add_argument('-p','--output_peptides', type = Path, help='Output Peptides')
    parser.add_argument('-r','--output_proteins', type = Path, help='Output Proteins')
    parser.add_argument('-c','--protein_coverage', type = Path, help='Added Protein Coverage')
    parser.add_argument('-d','--protein_coverage_external', type = Path, help='Added Protein Coverage (External)')
    parser.add_argument('-t','--cosine_cutoff', type = float, help='Cosine Cutoff')
    parser.add_argument('-l','--explained_intensity_cutoff', type = float, help='Explained Intensity Cutoff')
    parser.add_argument('-z','--nextprot_pe', type = Path, help='NextProt PEs (Latest)')
    parser.add_argument('-w','--nextprot_releases', type = Path, help='NextProt Releases')
    parser.add_argument('-m','--msv_to_pxd_mapping', type = Path, help='MSV to PXD Mapping')

    if len(sys.argv) < 4:
        parser.print_help()
        sys.exit(1)
    return parser.parse_args()

def add_brackets(pep):
    aa_breakpoints = []
    pep = pep.replace('+','[+').replace('-','[-')
    for i,aa in enumerate(pep[1:]):
        if not pep[i-1].isalpha() and pep[i].isalpha():
            aa_breakpoints.append(i)
        elif not pep[i].isalpha() and i == (len(pep)-2):
            aa_breakpoints.append(i+2)
    for i,breakpoint in enumerate(reversed(aa_breakpoints)):
        end_bracket = ']'
        if i == len(aa_breakpoints)-1 and pep[0] == '[':
            end_bracket = ']-'
        pep = pep[:breakpoint] + end_bracket + pep[breakpoint:]
    return pep

def msv_to_pxd(msv, msv_mapping):
    mapping = msv_mapping.get(msv,{}).get('px_accession')
    if not mapping:
        mapping = msv
    return mapping

def correct_usi(usi_input, msv_mapping):
    if '[' in usi_input or 'PXD' in usi_input:
        return usi_input
    else:
        split_usi = usi_input.split(':')
        split_usi[1] = msv_to_pxd(split_usi[1], msv_mapping)
        split_usi[5] = '/'.join([add_brackets(split_usi[5].split('/')[0]),split_usi[5].split('/')[1]])
    #     return usi_link(':'.join(split_usi))
        return ':'.join(split_usi)


def find_overlap(existing_peptides, new_peptides, protein_length, protein_pe, name, max_seq_output = 10):
    kb_pos = set([sorted(positions, key = lambda x: x[0])[0] for positions in existing_peptides.values()])
    added_pos = set([sorted(positions, key = lambda x: x[0])[0] for positions in new_peptides.values()])

    kb_peptides = set(["{} ({}-{})".format(peptide.replace('L','I'),sorted(positions, key = lambda x: x[0])[0][0],sorted(positions, key = lambda x: x[0])[0][1]) for peptide, positions in existing_peptides.items()])
    added_peptides = set(["{} ({}-{})".format(peptide.replace('L','I'),sorted(positions, key = lambda x: x[0])[0][0],sorted(positions, key = lambda x: x[0])[0][1]) for peptide, positions in new_peptides.items()])

    # print(kb_peptides)
    # print(added_peptides)

    novel_peptides = added_peptides.difference(kb_peptides)
    supporting_peptides = added_peptides.intersection(kb_peptides)

    nonoverlapping_kb_peptides = hupo_nonoverlapping(kb_pos.difference(added_pos))
    nonoverlapping_added_peptides = hupo_nonoverlapping(added_pos.difference(kb_pos))
    nonoverlapping_intersection = hupo_nonoverlapping(kb_pos.union(added_pos))
    nonoverlapping_all_kb = hupo_nonoverlapping(kb_pos)

    added_coverage = find_coverage(list(existing_peptides.values()), list(new_peptides.values()), protein_length)
    total_coverage = find_coverage([],list(existing_peptides.values()) + list(new_peptides.values()), protein_length)

    return {
        'supporting_peptides'+name:','.join(list(supporting_peptides)[:max_seq_output]) if len(supporting_peptides) >= 1 else ' ',
        'novel_peptides'+name:','.join(list(novel_peptides)[:max_seq_output]) if len(novel_peptides) >= 1 else ' ',
        'kb_hpp'+name:nonoverlapping_all_kb,
        'new_hpp'+name:nonoverlapping_added_peptides,
        'combined_hpp'+name:nonoverlapping_intersection,
        'added_kb_coverage'+name:added_coverage,
        'total_coverage'+name:total_coverage,
        'coverage_increase'+name:total_coverage/(total_coverage-added_coverage) if total_coverage-added_coverage != 0 else 1000,
        'promoted'+name:'Yes' if (nonoverlapping_all_kb < 2 and nonoverlapping_intersection >= 2) else 'No'
    }, [s.split(' ')[0] for s in supporting_peptides]

def find_coverage(current_kb, new_kb, length):

    current_coverage = [False for _ in range(length)]
    added_coverage = [False for _ in range(length)]

    for p in current_kb:
        for (start,end) in p:
            for i in range(start,end+1):
                if i < len(current_coverage)-1:
                    current_coverage[i-1] = True
    for p in new_kb:
        for (start,end) in p:
            for i in range(start,end+1):
                if i < len(added_coverage)-1:
                    added_coverage[i-1] = True if not current_coverage[i-1] else False

    return len([c for c in added_coverage if c])

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

def read_protein_coverage(protein_coverage_file,all_proteins,added_proteins,pep_info,nextprot_pe):
    with open(protein_coverage_file) as f:
        r = csv.DictReader(f, delimiter='\t')
        for l in r:
            proteins_w_coords = l['proteins_w_coords']
            peptide = l['peptide']
            il_peptide = l['il_peptide']
            gene_unique = l['gene_unique'] == 'True'
            sp_matches_saav = int(l['num_proteins_no_iso_no_tr'])

            proteins = [
                [p.split(' ')[0]] + p.split(' ')[1].replace('(','').replace(')','').split('-')
                for p in
                proteins_w_coords.split(';')
            ]

            match_obj = {
                'gene_unique': gene_unique,
                'canonical_matches': int(sp_matches_saav),
                'all_proteins_w_coords': proteins_w_coords
            }

            pep_info[il_peptide] = match_obj

            for protein in proteins:
                if len(protein) == 3:
                    accession = protein[0].split('|')[1]
                    aa_start = int(protein[1])
                    aa_end = int(protein[2])
                    if accession in nextprot_pe:
                        if 'XXX_' in protein[0]:
                            accession = 'XXX_{}'.format(accession)
                        all_proteins[accession][il_peptide].append((aa_start,aa_end))
                        if (aa_end - aa_start >= 9) and gene_unique and sp_matches_saav == 1:
                            added_proteins[accession][il_peptide].append((aa_start,aa_end))
    return all_proteins,added_proteins,pep_info

def main():
    args = arguments()

    representative_per_precursor = {}
    protein_length = {}

    # note that this has a limit of 100000000 so with multi-species in ProteinExplorer this will likely need to change
    url = "http://massive.ucsd.edu/ProteoSAFe/ProteinLibraryServlet?task=protein_explorer_proteins&file=&pageSize=100000000&offset=0&query=%2523%257B%2522unroll%2522%253A%2522no%2522%252C%2522include_synthetics%2522%253A%2522yes%2522%252C%2522datasets%2522%253A%2522%2522%252C%2522accession_input%2522%253A%2522%2522%257D&query_type=representative&_=1560375217897"
    proteins = requests.get(url).json()['row_data']
    for protein in proteins:
        # nextprot_pe[protein['accession']] = int(protein['proteinexistance'])
        protein_length[protein['accession']] = int(protein['length'])
    kb_proteins = defaultdict(lambda: defaultdict(list))
    added_proteins = defaultdict(lambda: defaultdict(list))
    all_proteins = defaultdict(lambda: defaultdict(list))
    pep_mapping_info = {}
    kb_proteins_w_synthetic = defaultdict(lambda: defaultdict(list))
    added_proteins_w_synthetic = defaultdict(lambda: defaultdict(list))
    has_synthetic = set()

    added_proteins_matching_synthetic = defaultdict(lambda: defaultdict(list))
    added_proteins_explained_intensity = defaultdict(lambda: defaultdict(list))

    with open(args.msv_to_pxd_mapping) as json_file:
        msv_mapping = json.load(json_file)

    kb_seq = set()
    kb_seq_matching = set()

    with open(args.kb_pep) as f:
        r = csv.DictReader(f, delimiter='\t')
        for l in r:
            if int(l['library']) == 2:
                kb_seq.add(l['demodified'].replace('I','L'))
                kb_proteins[l['protein']][l['demodified'].replace('I','L')].append((int(l['aa_start'])-1,int(l['aa_end'])))
            if int(l['library']) == 3 or int(l['library']) == 4:
                has_synthetic.add(l['demodified'].replace('I','L'))

    kb_seq_has_synthetic = kb_seq.intersection(has_synthetic)

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

    nextprot_releases_pe = {}

    for nextprot_release in args.nextprot_releases.glob('*'):
        release_str = nextprot_release.stem
        nextprot_releases_pe[release_str] = {}
        for pe_file in nextprot_release.glob('ac_lists/*PE*'):
            pe = pe_file.stem.split('_')[3].replace('PE','')
            with open(pe_file) as f:
                for l in f:
                    nextprot_releases_pe[release_str][l.rstrip().replace('NX_','')] = pe

    all_proteins,added_proteins,pep_mapping_info = read_protein_coverage(args.protein_coverage,all_proteins,added_proteins,pep_mapping_info,nextprot_pe)

    for protein_coverage_file in args.protein_coverage_external.glob('*'):
        all_proteins,added_proteins,pep_mapping_info = read_protein_coverage(protein_coverage_file,all_proteins,added_proteins,pep_mapping_info,nextprot_pe)

    for protein in all_proteins:
        for peptide in all_proteins[protein].keys():
            if peptide.replace('I','L') in peptide_to_protein:
                peptide_to_protein[peptide.replace('I','L')] += ', ' + protein
            else:
                peptide_to_protein[peptide.replace('I','L')] = protein

    for protein in added_proteins:
        for peptide in added_proteins[protein].keys():
            peptide_to_protein_added[peptide.replace('I','L')] = protein
            if peptide.replace('I','L') in has_synthetic:
                added_proteins_w_synthetic[protein][peptide] = added_proteins[protein][peptide]

    with open(args.output_psms,'w') as w:
        header = ['protein','decoy','pe','ms_evidence','filename','scan','sequence','charge','usi','score','pass','type','parent_mass','frag_tol','synthetic_filename','synthetic_scan','synthetic_usi','cosine','synthetic_match','explained_intensity','hpp_match','gene_unique','canonical_matches','all_proteins_w_coords','aa_start','aa_end']

        o = csv.DictWriter(w, delimiter='\t',fieldnames = header)
        o.writeheader()
        for input_psm in chain(args.input_psms.glob('*'),args.input_psms_external.glob('*')):
            with open(input_psm) as f:
                r = csv.DictReader(f, delimiter='\t')
                for l in r:
                    il_peptide = ''.join([p.replace('I','L') for p in l['sequence'] if p.isalpha()])
                    protein = peptide_to_protein.get(il_peptide)
                    l['usi'] = correct_usi(l['usi'], msv_mapping)
                    l['synthetic_usi'] = correct_usi(l['synthetic_usi'], msv_mapping) if l['synthetic_usi'] != '' else l['synthetic_usi']
                    # l['sequence'] = l['sequence'] if '[' in l['sequence'] else add_brackets(l['sequence'])
                    if protein:
                        protein_no_decoy = protein.replace('XXX_','')
                        l.update({
                            'protein': protein,
                            'pe': nextprot_pe.get(protein_no_decoy,''),
                            'ms_evidence':ms_existence.get(protein_no_decoy,'no')
                        })
                        l.update(pep_mapping_info[il_peptide])
                    else:
                        l.update({
                            'protein': '',
                            'pe': '',
                            'ms_evidence':''
                        })
                    if il_peptide in all_proteins.get(protein,{}):
                        l['aa_start'],l['aa_end'] = all_proteins[protein][il_peptide][0]
                    else:
                        l['aa_start'],l['aa_end'] = "N/A","N/A"

                    if il_peptide in added_proteins[protein]:
                        if il_peptide in kb_seq:
                            l['type'] = 'Matches existing evidence'
                        else:
                            l['type'] = 'New protein evidence'
                        l['hpp_match'] = 'True'
                    else:
                        l['type'] = 'N/A'
                        l['hpp_match'] = 'False'
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
        header = ['protein','pe','ms_evidence','aa_total','database_filename','database_scan','database_usi','sequence','charge','score','pass','type','parent_mass','cosine_filename','cosine_scan','cosine_usi','synthetic_filename','synthetic_scan','synthetic_usi','cosine','synthetic_match','cosine_score_match','explained_intensity','hpp_match','gene_unique','canonical_matches','all_proteins_w_coords','aa_start','aa_end','frag_tol']
        r = csv.DictWriter(w, delimiter = '\t', fieldnames = header)
        r.writeheader()
        for (sequence, charge), best_psm in representative_per_precursor.items():
            sequence_il = ''.join([p.replace('I','L') for p in sequence if p.isalpha()])
            protein = peptide_to_protein.get(sequence_il)
            if protein:
                protein_no_decoy = protein.replace('XXX_','')
                best_psm.update({
                    'protein': protein,
                    'pe': nextprot_pe.get(protein_no_decoy,''),
                    'ms_evidence':ms_existence.get(protein_no_decoy,'no'),
                    'aa_total':protein_length.get(protein_no_decoy,0)
                })
                best_psm.update(pep_mapping_info[sequence_il])
            else:
                best_psm.update({
                    'protein': '',
                    'pe': '',
                    'ms_evidence':'',
                    'aa_total':''
                })
            if sequence_il in all_proteins.get(protein,{}):
                best_psm['aa_start'],best_psm['aa_end'] = all_proteins[protein][sequence_il][0]
            else:
                best_psm['aa_start'],best_psm['aa_end'] = "N/A","N/A"
            if sequence_il in added_proteins.get(protein,{}):
                if sequence_il in kb_seq:
                    best_psm['type'] = 'Matches existing evidence'
                else:
                    best_psm['type'] = 'New protein evidence'
                best_psm['hpp_match'] = 'True'
            else:
                best_psm['hpp_match'] = 'False'
                l['type'] = 'N/A'
            r.writerow(best_psm)
            if sequence_il in added_proteins[protein]:
                if best_psm['explained_intensity'] >= args.explained_intensity_cutoff:
                    if best_psm['cosine'] >= args.cosine_cutoff:
                        added_proteins_matching_synthetic[protein][sequence_il] = added_proteins[protein][sequence_il]
                    added_proteins_explained_intensity[protein][sequence_il] = added_proteins[protein][sequence_il]

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
            'total_coverage',
            'coverage_increase',
            'supporting_peptides_w_synthetic',
            'novel_peptides_w_synthetic',
            'kb_hpp_w_synthetic',
            'new_hpp_w_synthetic',
            'combined_hpp_w_synthetic',
            'added_kb_coverage_w_synthetic',
            'total_coverage_w_synthetic',
            'coverage_increase_w_synthetic',
            'supporting_peptides_without_kb',
            'novel_peptides_without_kb',
            'kb_hpp_without_kb',
            'new_hpp_without_kb',
            'combined_hpp_without_kb',
            'added_kb_coverage_without_kb',
            'total_coverage_without_kb',
            'coverage_increase_without_kb',
            'promoted',
            'promoted_w_synthetic',
            'promoted_without_kb',
            'cosine_cutoff',
            'explained_intensity_cutoff'
        ] + ['_dyn_#{}'.format(release) for release in sorted(nextprot_releases_pe.keys())]

        w = csv.DictWriter(fo, delimiter = '\t', fieldnames = fieldnames)
        w.writeheader()

        for protein in set(all_proteins.keys()).union(set(kb_proteins.keys())):
            protein_no_decoy = protein.replace('XXX_','')
            if (protein_no_decoy in nextprot_pe):
                protein_dict = {
                    'protein': protein,
                    'pe': nextprot_pe.get(protein_no_decoy,0),
                    'ms_evidence':ms_existence.get(protein_no_decoy,'no'),
                    'aa_total':protein_length.get(protein_no_decoy,0)
                }
                protein_dict.update(find_overlap(kb_proteins[protein],added_proteins_explained_intensity[protein],int(protein_dict['aa_total']),int(protein_dict['pe']),'')[0])
                protein_dict.update(find_overlap(kb_proteins_w_synthetic[protein],added_proteins_matching_synthetic[protein],int(protein_dict['aa_total']),int(protein_dict['pe']),'_w_synthetic')[0])
                protein_dict.update(find_overlap({},added_proteins_explained_intensity[protein],int(protein_dict['aa_total']),int(protein_dict['pe']),'_without_kb')[0])
                protein_dict.update({
                    'cosine_cutoff':args.cosine_cutoff,
                    'explained_intensity_cutoff':args.explained_intensity_cutoff
                })
                for release, pe_dict in nextprot_releases_pe.items():
                    protein_dict['_dyn_#{}'.format(release)] = pe_dict.get(protein_no_decoy,0)
                w.writerow(protein_dict)

if __name__ == '__main__':
    main()
