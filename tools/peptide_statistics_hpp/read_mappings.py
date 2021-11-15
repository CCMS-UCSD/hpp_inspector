import csv
import sys
from collections import defaultdict
from python_ms_utilities import mapping
import pandas as pd
from datetime import datetime

csv.field_size_limit(sys.maxsize)

def read_protein_coverage(protein_coverage_file,seen_sequences,proteome, filter = None, output_representatives = True, use_job_level_thresholds = False):

    pep_info = {}
    protein_mapping_out = defaultdict(lambda: defaultdict(set))
    peptide_to_exon_map = defaultdict(list)
    output_peptides = []

    already_seen = 0
    need_to_parse = 0

    protein_hpp_fdr = {}
    protein_hint_fdr = {}

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

                all_protein_fdr = float(l.get('hint_protein_fdr',-1))
                hpp_protein_fdr = float(l.get('hpp_protein_fdr',-1))
                cosine = float(l.get('cosine',-1))

                pass_filters = True

                if filter:
                    pass_precursor_fdr = float(l.get('precursor_fdr',-1)) <= float(l.get('precursor_fdr_cutoff',-1))
                    pass_cosine = cosine >= float(l.get('cosine_cutoff',-1))
                    pass_expl_intensity = float(l.get('explained_intensity',1)) >= float(l.get('explained_intensity_cutoff',0))
                    pass_annotated_ions = float(l.get('matched_ions',1000)) >= float(l.get('matched_ions_cutoff',0))
                    # print(pass_precursor_fdr,pass_cosine,pass_expl_intensity,pass_annotated_ions)
                    pass_filters = pass_precursor_fdr and pass_annotated_ions and (pass_cosine or pass_expl_intensity)

                if use_job_level_thresholds:
                    cosine = 1 if cosine >= float(l.get('cosine_cutoff',-1)) else cosine

                protein_mappings = mapping.string_to_protein_mappings(mapped_protein_str)
                exon_mappings = mapping.string_to_exon_mappings(mapped_exon_str)

                m = mapping.summarize_protein_mappings(proteome,protein_mappings,exon_mappings)

                gene_unique = m['gene_unique_incl_mismatch']
                sp_matches_saav = int(m['num_proteins_no_iso_no_tr_incl_mismatch'])

                match_obj = {
                    'gene_unique': gene_unique,
                    'canonical_matches': int(sp_matches_saav),
                    'len': len(il_peptide),
                    # remove gene-uniqueness from HPP criteria
                    'hpp': len(il_peptide) >= 9 and sp_matches_saav == 1,
                    'all_proteins_w_coords': mapped_protein_str,
                    'all_mapped_exons': mapped_exon_str,
                    'exon_junctions_covered': m['exon_junctions_covered'],
                    'exons_covered_no_junction': m['exons_covered_no_junction'],
                    'total_unique_exons_covered': m['total_unique_exons_covered'],
                    'mapped_proteins' : []
                }

                if pass_filters:

                    for protein_mapping in protein_mappings:
                        if len(protein_mapping.mismatches) == 0:
                            protein_mapping_out[protein_mapping.protein_accession][il_peptide].add((protein_mapping.start_pos,protein_mapping.end_pos,cosine,all_protein_fdr,hpp_protein_fdr,len(il_peptide) >= 9 and sp_matches_saav == 1))
                            match_obj['mapped_proteins'].append((protein_mapping.protein_accession,(protein_mapping.start_pos,protein_mapping.end_pos),protein_mapping.il_ambiguous))
                            protein_hpp_fdr[protein_mapping.protein_accession] = min(protein_hpp_fdr.get(protein_mapping.protein_accession,1),hpp_protein_fdr)
                            protein_hint_fdr[protein_mapping.protein_accession] = min(protein_hint_fdr.get(protein_mapping.protein_accession,1),all_protein_fdr)
                    for exon_mapping in exon_mappings:
                        for (complete, mapped) in zip(exon_mapping.complete_coordinates, exon_mapping.matched_coordinates):
                            peptide_to_exon_map[(chr,complete,exon_mapping.gene,','.join(exon_mapping.transcripts))].append((peptide,mapped,exon_mapping.matched_coordinates,exon_mapping.complete_coordinates))

                    pep_info[peptide] = match_obj
            else:
                already_seen += 1
                
            if output_representatives:
                output_peptides.append(l)
            
        print("Seen {}, need to parse {} mappings".format(already_seen,need_to_parse))

    return protein_mapping_out,pep_info,peptide_to_exon_map,output_peptides,protein_hpp_fdr,protein_hint_fdr
