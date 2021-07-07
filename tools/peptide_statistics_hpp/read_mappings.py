import csv
import sys
from collections import defaultdict
from python_ms_utilities import mapping
import pandas as pd
from datetime import datetime

csv.field_size_limit(sys.maxsize)

def read_protein_coverage(protein_coverage_file,seen_sequences,proteome, filter = False, output_representatives = True):

    pep_info = {}
    protein_mapping_out = defaultdict(lambda: defaultdict(set))
    peptide_to_exon_map = defaultdict(list)
    output_peptides = []

    already_seen = 0
    need_to_parse = 0

    current_time = datetime.now()

    with open(protein_coverage_file) as f:
        r = csv.DictReader(f, delimiter='\t')

        for l in r:

            peptide = l['peptide'] if 'peptide' in l else l['sequence_unmodified']

            if peptide not in seen_sequences:
                need_to_parse += 1
                il_peptide = l['il_peptide'] if 'il_peptide' in l else l['sequence_unmodified_il']

                mapped_protein_str = l['mapped_proteins'] if 'mapped_proteins' in l else l['all_proteins_w_coords']
                mapped_exon_str = l['mapped_exons'] if 'mapped_exons' in l else l['all_mapped_exons']

                precursor_fdr = float(l.get('precursor_fdr',-1))
                expl_intensity = float(l.get('explained_intensity',1))
                annoated_ions = float(l.get('matched_ions',1000))

                protein_mappings = mapping.string_to_protein_mappings(mapped_protein_str)
                exon_mappings = mapping.string_to_exon_mappings(mapped_exon_str)

                m = mapping.summarize_protein_mappings(proteome,protein_mappings,exon_mappings)

                gene_unique = m['gene_unique_incl_mismatch']
                sp_matches_saav = int(m['num_proteins_no_tr_incl_mismatch'])

                match_obj = {
                    'gene_unique': gene_unique,
                    'canonical_matches': int(sp_matches_saav),
                    'len': len(il_peptide),
                    'hpp': len(il_peptide) >= 9 and gene_unique and sp_matches_saav == 1,
                    'all_proteins_w_coords': mapped_protein_str,
                    'all_mapped_exons': mapped_exon_str,
                    'exon_junctions_covered': m['exon_junctions_covered'],
                    'exons_covered_no_junction': m['exons_covered_no_junction'],
                    'total_unique_exons_covered': m['total_unique_exons_covered'],
                    'mapped_proteins' : []
                }

                if not filter or (precursor_fdr <= 0.01 and expl_intensity >= 0.4 and annoated_ions >= 6):

                    for protein_mapping in protein_mappings:
                        if len(protein_mapping.mismatches) == 0:
                            protein_mapping_out[protein_mapping.protein_accession][il_peptide].add((protein_mapping.start_pos,protein_mapping.end_pos,l.get('cosine')))
                            match_obj['mapped_proteins'].append((protein_mapping.protein_accession,(protein_mapping.start_pos,protein_mapping.end_pos),protein_mapping.il_ambiguous))
                    for exon_mapping in exon_mappings:
                        for (complete, mapped) in zip(exon_mapping.complete_coordinates, exon_mapping.matched_coordinates):
                            peptide_to_exon_map[(chr,complete,exon_mapping.gene,','.join(exon_mapping.transcripts))].append((peptide,mapped,exon_mapping.matched_coordinates,exon_mapping.complete_coordinates))

                    pep_info[peptide] = match_obj
            else:
                already_seen += 1
                
            if output_representatives:
                output_peptides.append(l)
            
        print("Seen {}, need to parse {} mappings".format(already_seen,need_to_parse))

    return protein_mapping_out,pep_info,peptide_to_exon_map,output_peptides
