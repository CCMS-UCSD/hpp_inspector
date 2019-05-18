import argparse
import sys
import csv
import glob
import mztab


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

def arguments():
    parser = argparse.ArgumentParser(description='mzTab to list of peptides')
    parser.add_argument('-p','--kb_pep', type = str, help='Peptides from KB')
    parser.add_argument('-c','--protein_coverage', type = str, help='Added Protein Coverage')
    parser.add_argument('-n','--novel_coverage', type = str, help='Novel Coverage')
    parser.add_argument('-x','--nextprot_pe', type = str, help='NextProt PEs')

    if len(sys.argv) < 3:
        parser.print_help()
        sys.exit(1)
    return parser.parse_args()

def main():
    args = arguments()
    nextprot_pe = defaultdict(lambda: 0)
    with open(args.nextprot_pe) as f:
        r = csv.DictReader(f, delimiter='\t')
        for l in r:
            nextprot_pe[l['Protein'].replace('NX_','')] = int(l['PE'])
    kb_proteins = defaultdict(lambda: defaultdict(list))
    added_proteins = defaultdict(lambda: defaultdict(list))
    with open(args.kb_pep) as f:
        r = csv.DictReader(f, delimiter='\t')
        for l in r:
            kb_proteins[l['protein']][l['peptide']].append((int(l['aa_start']),int(l['aa_end'])))
    with open(args.protein_coverage) as f:
        r = csv.DictReader(f, delimiter='\t')
        for l in r:
            proteins_w_coords = l['proteins_w_coords']
            peptide = l['peptide']
            gene_unique = l['gene_unique'] == 'True'
            if gene_unique:
                proteins = [
                    [p.split(' ')[0]] + p.split(' ')[1].replace('(','').replace(')','').split('-')
                    for p in
                    proteins_w_coords.split(';')
                ]
                for protein in proteins:
                    accession = protein[0].split('|')[1]
                    aa_start = int(protein[1])
                    aa_end = int(protein[2])
                    if nextprot_pe[accession] > 1:
                        added_proteins[accession][peptide].append((aa_start,aa_end))

    for protein, added_peptides in added_proteins.item():
        kb_peptides = kb_proteins[protein]
        hupo_nonoverlapping
        unique_kb_peptides = len(set(kb_peptides.keys()).difference(set(added_peptides.keys())))
        unique_added_peptides = len(set(added_peptides.keys()).difference(set(kb_peptides.keys())))
        unique_intersection = len(set(added_peptides.keys()).difference(set(kb_peptides.keys())))

        nonoverlapping_kb_peptides = len(set(kb_peptides.keys()).difference(set(added_peptides.keys())))
        nonoverlapping_added_peptides = len(set(added_peptides.keys()).difference(set(kb_peptides.keys())))
        nonoverlapping_intersection = len(set(added_peptides.keys()).difference(set(kb_peptides.keys())))

if __name__ == '__main__':
    main()
