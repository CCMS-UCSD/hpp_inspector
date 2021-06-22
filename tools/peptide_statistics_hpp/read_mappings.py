import csv
import sys
from python_ms_utilities import mapping

csv.field_size_limit(sys.maxsize)

def read_protein_coverage(protein_coverage_file,all_proteins,added_proteins,pep_info,peptide_to_exon_map,proteome):
    with open(protein_coverage_file) as f:
        r = csv.DictReader(f, delimiter='\t')
        for l in r:
            peptide = l['peptide']
            il_peptide = l['il_peptide']

            protein_mappings = mapping.string_to_protein_mappings(l['mapped_proteins'])
            exon_mappings = mapping.string_to_exon_mappings(l['mapped_exons'])

            m = mapping.summarize_protein_mappings(proteome,protein_mappings,exon_mappings)

            gene_unique = m['gene_unique_incl_mismatch']
            sp_matches_saav = int(m['num_proteins_no_tr_incl_mismatch'])

            match_obj = {
                'gene_unique': gene_unique,
                'canonical_matches': int(sp_matches_saav),
                'all_proteins_w_coords': l['mapped_proteins'],
                'all_mapped_exons': l['mapped_exons'],
                'exon_junctions_covered': m['exon_junctions_covered'],
                'exons_covered_no_junction': m['exons_covered_no_junction'],
                'total_unique_exons_covered': m['total_unique_exons_covered'],
                'mapped_proteins' : []
            }

            for protein_mapping in protein_mappings:
                if len(protein_mapping.mismatches) == 0:
                    match_obj['mapped_proteins'].append((protein_mapping.protein_accession,(protein_mapping.start_pos,protein_mapping.end_pos),protein_mapping.il_ambiguous))
                    all_proteins[protein_mapping.protein_accession][il_peptide].append(((protein_mapping.start_pos,protein_mapping.end_pos),gene_unique and sp_matches_saav == 1))
                    if (protein_mapping.end_pos - protein_mapping.start_pos) >= 8 and gene_unique and sp_matches_saav == 1:
                        added_proteins[protein_mapping.protein_accession][il_peptide].append((protein_mapping.start_pos,protein_mapping.end_pos))

            for exon_mapping in exon_mappings:
                for (complete, mapped) in zip(exon_mapping.complete_coordinates, exon_mapping.matched_coordinates):
                    peptide_to_exon_map[(chr,complete,exon_mapping.gene,','.join(exon_mapping.transcripts))].append((peptide,mapped,exon_mapping.matched_coordinates,exon_mapping.complete_coordinates))

            pep_info[peptide] = match_obj

    return all_proteins,added_proteins,pep_info,peptide_to_exon_map
