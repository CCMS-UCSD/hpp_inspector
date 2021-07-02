import sys
import csv
from csv import DictReader, DictWriter
from collections import defaultdict
from python_ms_utilities import fdr

csv.field_size_limit(sys.maxsize)

def non_nested_score(segments):
    if len(segments) == 0:
        return 0
    max_peps = [[0 for _ in range(2)] for _ in range(len(segments))]
    max_score_at_pep = [[0 for _ in range(2)] for _ in range(len(segments))]
    max_peps_count = [[0 for _ in range(2)] for _ in range(len(segments))]
    segments = [(int(s[0]),int(s[1]),float(s[2])) for s in segments]
    ordered_segments = sorted(segments, key = lambda x: (int(x[1]),float(x[2]),int(x[0])))
    for i in range(len(segments)):
        max_peps_j = 0
        max_score_at_pep_j = 0
        j = -1
        for k in range(i-1,-1,-1):
            if (ordered_segments[k][0] < ordered_segments[i][0] and
                ordered_segments[k][1] < ordered_segments[i][1] and
                ordered_segments[k][1] >= ordered_segments[i][0]):
                if j >= 0:
                    if (max(max_peps[k]) > max(max_peps[j])):
                        j = k
                else:
                    j = k
        if j >= 0:
            max_peps_j = max(max_peps[j])
        j = -1
        for k in range(i-1,-1,-1):
            if (ordered_segments[k][1] <= ordered_segments[i][0]):
                j = k
                break
        if j >= 0:
            max_score_at_pep_j = max(max_score_at_pep[j])
        if (max(max_peps_j,max_score_at_pep_j)) == 0:
            max_peps[i][0] = ordered_segments[i][2]
        else:
            max_peps[i][1] = ordered_segments[i][2] + max(max_peps_j,max_score_at_pep_j)
        max_score_at_pep_j = 0
        j = -1
        for k in range(i-1,-1,-1):
            if (ordered_segments[k][1] < ordered_segments[i][1] and
                ordered_segments[k][1] >= ordered_segments[i][0]):
                if j >= 0:
                    if (max(max_score_at_pep[k]) > max(max_score_at_pep[j])):
                        j = k
                else:
                    j = k
        max_score_at_pep_j_0 = 0
        max_score_at_pep_j_1 = 0
        if j >= 0:
            max_score_at_pep_j_0 = max_score_at_pep[j][0]
            max_score_at_pep_j_1 = max_score_at_pep[j][1]
        max_score_at_pep_i_prev_0 = 0
        max_score_at_pep_i_prev_1 = 0
        max_peps_count_prev_i = 0
        if i != 0:
            max_score_at_pep_i_prev_0 = max_peps[i-1][0]
            max_score_at_pep_i_prev_1 = max_peps[i-1][1]
            max_peps_count_prev_i = max_peps_count[i-1][1]
        max_score_at_pep[i][0] = max(max_score_at_pep_j_0,max_peps[i][0],max_score_at_pep_i_prev_0)
        max_score_at_pep[i][1] = max(max_score_at_pep_j_1,max_peps[i][1],max_score_at_pep_i_prev_1)
        max_peps_count[i][1] = 1+max_peps_count_prev_i if max_score_at_pep[i][1] == max_peps[i][1] and max_peps[i][1] != max_score_at_pep_i_prev_1 else max_peps_count_prev_i
    return max_score_at_pep[-1][-1], max_peps_count[-1][-1]

input = sys.argv[1]
output = sys.argv[2]

precursors_per_protein = defaultdict(lambda: defaultdict(float))

with open(input) as f:
    r = DictReader(f, delimiter = '\t')
    for l in r:
        proteins = l['protein'].split(' ###')[0].split(';')
        if 'Canonical' in l['protein_type'] and len(proteins) == 1 and l['hpp_match'] == 'True':
            pos = (int(l['aa_start']),int(l['aa_end']))
            sequence = ''.join([s.replace('I','L') for s in l['sequence'] if s.isalpha()])
            prev_score = precursors_per_protein[proteins[0]][pos]
            precursors_per_protein[proteins[0]][pos] = max(float(l['score']),prev_score)

protein_w_scores = []

for protein, precursors in precursors_per_protein.items():
    score, count = mapping.non_nested_score([(*k,v) for k,v in precursors.items()])
    protein_w_scores.append(fdr.ScoredElement(protein,'XXX_' in protein, score))

fdr_dict = fdr.calculate_fdr(protein_w_scores)

with open(output, 'w') as f:
    w = DictWriter(f, delimiter = '\t', fieldnames = ['Protein','Score','NumPeptides','FDR'])
    w.writeheader()
    for protein,decoy,score in protein_w_scores:
        w.writerow({
            'Protein':protein,
            'Score':score,
            'NumPeptides':0,
            'FDR':fdr_dict[protein]
        })
