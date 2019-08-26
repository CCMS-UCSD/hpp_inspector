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
import pathlib
import pickle


CUTOFF = {"6": 2.12167006027, "7": 2.32154253572, "8": 2.53031824698, "9": 2.7717776307, "10": 3.0419577266, "11": 3.20404758311, "12": 3.4028353721, "13": 3.14854063781, "14": 3.47053358426, "15": 3.35416414915, "16": 3.34418455532, "17": 3.24904542906, "18": 3.18426919317, "19": 3.01362943461, "20": 3.12316407632, "21": 3.08160158158, "22": 2.59466460406, "23": 2.96230440256, "24": 3.11610830789, "25": 2.93990420679, "26": 2.6947192629, "27": 2.43042157531, "28": 2.5287628139300002, "29": 2.26034401643, "30": 2.60979068254, "31": 2.70004841417, "32": 2.69919517925, "33": 2.18110553715, "34": 1.90115111418, "35": 1.5402648112000001, "36": 1.74919562203, "37": 1.88066887473, "38": 1.58471929702, "39": 1.73377627878, "40": 3.10312899149}

def arguments():
    parser = argparse.ArgumentParser(description='mzTab to list of peptides')
    parser.add_argument('-p','--kb_pep', type = str, help='Peptides from KB')
    parser.add_argument('-m','--mztab', type = str, help='mzTab')
    parser.add_argument('-c','--protein_coverage', type = str, help='Added Protein Coverage')
    parser.add_argument('-s','--spec_on_server', type = str, help='Spectra')
    parser.add_argument('-n','--novel_proteins', type = str, help='Novel Proteins')
    parser.add_argument('-x','--novel_psms', type = str, help='Novel PSMs')
    parser.add_argument('-e','--novel_peptides', type = str, help='Novel Peptides')
    parser.add_argument('-y','--synthetics', type = str, help='Matched Synthetics')
    parser.add_argument('-t','--peak_tolerance', type = str, help='Peak Tolerance for Matching')
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
    nextprot_pe = defaultdict(lambda: 0)
    protein_length = defaultdict(lambda: 0)
    with open(args.synthetics,'rb') as f:
        synthetic_scans = pickle.load(f)

    tol = float(args.peak_tolerance)

    ids = {}

    if args.mztab != "" and args.mztab != None:
        for mztab_file in glob.glob(args.mztab + '/*'):
            try:
                ids = mztab.read(mztab_file, ids)
            except:
                ids = mztab.read_lib(mztab_file, ids)

    peptide_to_psm = defaultdict(list)
    for filescan,psm in ids.items():
        peptide_to_psm[''.join([p.replace('I','L') for p in psm[0].sequence if p.isalpha()])].append((filescan,psm))

    # print(peptide_to_psm)

    supporting_peptides = set()

    peptide_to_protein = {}

    psms_to_consider = defaultdict(set)

    ms_existence = {}

    cosine_to_synthetic = defaultdict(lambda: (0,('N/A','N/A')))

    with open(args.nextprot_pe) as f:
        r = csv.DictReader(f, delimiter='\t')
        for l in r:
            ms_existence[l['protein'].replace('NX_','')] = l['ms_evidence']

    # note that this has a limit of 100000000 so with multi-species in ProteinExplorer this will likely need to change
    url = "http://ccms-internal.ucsd.edu/ProteoSAFe/ProteinLibraryServlet?task=protein_explorer_proteins&file=&pageSize=100000000&offset=0&query=%2523%257B%2522unroll%2522%253A%2522no%2522%252C%2522include_synthetics%2522%253A%2522yes%2522%252C%2522datasets%2522%253A%2522%2522%252C%2522accession_input%2522%253A%2522%2522%257D&query_type=representative&_=1560375217897"
    proteins = requests.get(url).json()['row_data']
    for protein in proteins:
        nextprot_pe[protein['accession']] = int(protein['proteinexistance'])
        protein_length[protein['accession']] = int(protein['length'])
    kb_proteins = defaultdict(lambda: defaultdict(list))
    added_proteins = defaultdict(lambda: defaultdict(list))
    kb_proteins_w_synthetic = defaultdict(lambda: defaultdict(list))
    added_proteins_w_synthetic = defaultdict(lambda: defaultdict(list))
    has_synthetic = set()

    with open(args.kb_pep) as f:
        r = csv.DictReader(f, delimiter='\t')
        for l in r:
            if int(l['library']) == 2:
                kb_proteins[l['protein']][l['demodified']].append((int(l['aa_start']),int(l['aa_end'])))
            if int(l['library']) == 3 or int(l['library']) == 4:
                has_synthetic.add(l['demodified'])
        for protein in kb_proteins:
            for seq in kb_proteins[protein].keys():
                if seq in has_synthetic:
                    kb_proteins_w_synthetic[protein][seq] = kb_proteins[protein][seq]


    with open(args.protein_coverage) as f:
        r = csv.DictReader(f, delimiter='\t')
        for l in r:
            proteins_w_coords = l['proteins_w_coords']
            peptide = l['il_peptide']
            gene_unique = l['gene_unique'] == 'True'

            if gene_unique:
                proteins = [
                    [p.split(' ')[0]] + p.split(' ')[1].replace('(','').replace(')','').split('-')
                    for p in
                    proteins_w_coords.split(';')
                ]
                for protein in proteins:
                    if len(protein) == 3:
                        accession = protein[0].split('|')[1]
                        aa_start = int(protein[1]) + 1
                        aa_end = int(protein[2])
                        if accession in nextprot_pe and (aa_end - aa_start >= 9):
                            added_proteins[accession][peptide].append((aa_start,aa_end))

    for protein in added_proteins:
        for peptide in added_proteins[protein].keys():
            peptide_to_protein[peptide.replace('I','L')] = protein
            for (filescan,psm) in peptide_to_psm[peptide.replace('I','L')]:
                if synthetic_scans.get((psm[0].sequence,psm[0].charge)):
                    psms_to_consider[filescan[0]].add(filescan[1])
                    added_proteins_w_synthetic[protein][peptide] = added_proteins[protein][peptide]

    for filename in psms_to_consider:
        try:
            exts = pathlib.Path(filename).suffixes
            if 'MSV' in filename:
                filepath = filename.replace('f.','/data/massive/')
                if filepath[0] == 'M':
                    filepath = '/data/massive/' + filepath
            if exts[0] == '.mzML':
                with mzml.read(filepath) as reader:
                    for spectrum in reader:
                        peaks = zip(spectrum['m/z array'],spectrum['intensity array'])
                        spec_id = {idx.split('=')[0]:idx.split('=')[1] for idx in spectrum['id'].split(' ')}
                        scan = int(spec_id['scan'])
                        if scan in psms_to_consider[filename]:
                            synthetic_peaks,synthetic_filescan = synthetic_scans.get((ids[(filename,scan)][0].sequence,ids[(filename,scan)][0].charge))
                            cosine,_ = sa.score_alignment(list(peaks),list(synthetic_peaks),0,0,tol)
                            cosine_to_synthetic[(filename,scan)] = (cosine,synthetic_filescan)
        except:
            pass
        if exts[0] == '.mgf':
            pass
        if exts[0] == '.mzXML':
            pass


    with open(args.novel_proteins, 'w') as w:
        fieldnames = ['protein','gene','pe','aa_total','ms_evidence','fdr','supporting_peptides','novel_peptides','kb_hpp','new_hpp','combined_hpp','added_kb_coverage','supporting_peptides_w_synthetic','novel_peptides_w_synthetic','kb_hpp_w_synthetic','new_hpp_w_synthetic','combined_hpp_w_synthetic','added_kb_coverage_w_synthetic','promoted','promoted_w_synthetic']
        r = csv.DictWriter(w, delimiter = '\t', fieldnames = fieldnames)
        r.writeheader()

        for protein in added_proteins.keys():
            if (protein in nextprot_pe):
                protein_dict = {
                    'protein': protein,
                    'pe': nextprot_pe[protein],
                    'ms_evidence':ms_existence.get(protein,'no'),
                    'aa_total':protein_length.get(protein,0)
                }
                all_proteins_dict, all_supporting = find_overlap(kb_proteins[protein],added_proteins[protein],protein_length[protein],nextprot_pe[protein],'')
                protein_dict.update(all_proteins_dict)
                protein_dict.update(find_overlap(kb_proteins_w_synthetic[protein],added_proteins_w_synthetic[protein],protein_length[protein],nextprot_pe[protein],'_w_synthetic')[0])

                for supporting in all_supporting:
                    supporting_peptides.add(supporting)

                r.writerow(protein_dict)


    with open(args.novel_psms, 'w') as fw_psm, open(args.novel_peptides, 'w') as fw_pep:
        header = ['protein','pe','ms_evidence','filename','scan','sequence','charge','score','pass','type','parent_mass','synthetic_filename','synthetic_scan','best_cosine_to_synthetic']
        w_psm = csv.DictWriter(fw_psm, delimiter = '\t', fieldnames = header)
        w_pep = csv.DictWriter(fw_pep, delimiter = '\t', fieldnames = header)
        w_psm.writeheader()
        w_pep.writeheader()
        for peptide, protein in peptide_to_protein.items():
            # print(peptide,protein)
            best_psm = None
            best_psm_score = -1e9
            # print(peptide_to_psm[peptide.replace('I','L')])
            for (filescan,psm) in peptide_to_psm[peptide.replace('I','L')]:
                cosine, best_synthetic = cosine_to_synthetic[filescan]
                if best_synthetic[0] != 'N/A':
                    synthetic_filename = 'f.' + best_synthetic[0]
                else:
                    synthetic_filename = best_synthetic[0]
                psm_row = {
                    'protein':protein,
                    'ms_evidence':ms_existence.get(protein,'no'),
                    'pe':nextprot_pe.get(protein,0),
                    'filename':filescan[0],
                    'scan':filescan[1],
                    'sequence':psm[0].sequence,
                    'charge':psm[0].charge,
                    'score':psm[0].score,
                    'pass':'Above' if psm[0].score > CUTOFF[str(len(peptide))] else 'Below',
                    'parent_mass':psm[0].parent_mass,
                    'type':'Supporting' if peptide in supporting_peptides else 'Novel',
                    'synthetic_filename': synthetic_filename,
                    'synthetic_scan': best_synthetic[1],
                    'best_cosine_to_synthetic': cosine
                }
                w_psm.writerow(psm_row)
                if psm[0].score >= best_psm_score:
                    best_psm = psm_row
            if best_psm:
                w_pep.writerow(best_psm)

if __name__ == '__main__':
    main()
