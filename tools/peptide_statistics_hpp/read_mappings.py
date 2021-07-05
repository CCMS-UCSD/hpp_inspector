import csv
import sys
from python_ms_utilities import mapping

csv.field_size_limit(sys.maxsize)

def read_protein_coverage(protein_coverage_file,proteome):
    with open(protein_coverage_file) as f:
        r = csv.DictReader(f, delimiter='\t')

        peptide_info_out = []
        protein_mappings_out = []
        exon_mappings_out = []
        representatives_out = []

        for l in r:
            peptide = l['peptide'] if 'peptide' in l else l['sequence_unmodified']
            il_peptide = l['il_peptide'] if 'il_peptide' in l else l['sequence_unmodified_il']

            mapped_protein_str = l['mapped_proteins'] if 'mapped_proteins' in l else l['all_proteins_w_coords']
            mapped_exon_str = l['mapped_exons'] if 'mapped_exons' in l else l['all_mapped_exons']

            precursor_fdr = float(l.get('precursor_fdr',-1))

            protein_mappings = mapping.string_to_protein_mappings(mapped_protein_str)
            exon_mappings = mapping.string_to_exon_mappings(mapped_exon_str)

            m = mapping.summarize_protein_mappings(proteome,protein_mappings,exon_mappings)

            gene_unique = m['gene_unique_incl_mismatch']
            sp_matches_saav = int(m['num_proteins_no_tr_incl_mismatch'])

            match_obj = {
                'peptide': peptide,
                'gene_unique': gene_unique,
                'canonical_matches': int(sp_matches_saav),
                'all_proteins_w_coords': mapped_protein_str,
                'all_mapped_exons': mapped_exon_str,
                'exon_junctions_covered': m['exon_junctions_covered'],
                'exons_covered_no_junction': m['exons_covered_no_junction'],
                'total_unique_exons_covered': m['total_unique_exons_covered'],
            }

            if precursor_fdr <= 0.01:

                for protein_mapping in protein_mappings:
                    mapping_obj = {
                        'peptide':peptide,
                        'protein_accession':protein_mapping.protein_accession,
                        'start_pos':protein_mapping.start_pos,
                        'end_pos':protein_mapping.end_pos,
                        'il_ambiguous':protein_mapping.il_ambiguous,
                        'mismatches':protein_mapping.mismatches,
                        'synthetic_cosine':l.get('cosine'),
                    }
                    protein_mappings_out.append(mapping_obj)

                for exon_mapping in exon_mappings:
                    for (complete, mapped) in zip(exon_mapping.complete_coordinates, exon_mapping.matched_coordinates):
                        exon_mapping_obj = {
                            'chr':chr,
                            'complete_exon':complete,
                            'gene':exon_mapping.gene,
                            'all_transcripts':','.join(exon_mapping.transcripts),
                            'peptide':peptide,
                            'mapped':mapped,
                            'matched_coordinates':exon_mapping.matched_coordinates,
                            'complete_coordinates':exon_mapping.complete_coordinates
                        }
                        exon_mappings_out.append(exon_mapping_obj)

                peptide_info_out.append(match_obj)
                if 'database_filename' in l:
                    representatives_out.append({
                        'database_filename':l['database_filename'],
                        'database_scan':l['database_scan'],
                        'database_usi':l['database_usi'],
                        'cosine_filename':l['cosine_filename'],
                        'cosine_scan':l['cosine_scan'],
                        'cosine_usi':l['cosine_usi'],
                        'sequence':l['sequence'],
                        'charge':l['charge'],
                        'score':float(l['score']),
                        'parent_mass':l['parent_mass'],
                        'frag_tol':l['frag_tol'],
                        'synthetic_filename':l['synthetic_filename'],
                        'synthetic_scan':l['synthetic_scan'],
                        'synthetic_usi':l['synthetic_usi'],
                        'cosine':float(l['cosine']),
                        'synthetic_match':l['synthetic_match'],
                        'explained_intensity':float(l['explained_intensity']),
                        'matched_ions':l['matched_ions']
                    })

    return peptide_info_out,protein_mappings_out,exon_mappings_out,representatives_out
