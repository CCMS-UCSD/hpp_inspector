
def headers(key):
    return {
        'proteins':
            [
             'accession',
             'gene',
             'description',
             'uniprot_pe',
             'nextprot_pe',
             'sequence',
             'sequence_il',
             'annotations_rendered',
             'version',
             'species'
            ],
        'variants':
            [
             'id',
             'sequence',
             'demodified',
             'peptide_id',
             'precursor',
             'modifications',
            ],
        'peptides':
            [
             'id',
             'demodified',
             'number_of_proteins',
             'number_of_proteins_mismatch',
             'exon_unique',
             'splice_junction',
             'exon_mapped',
             'length'
            ],
        'peptide_mapping':
            [
             'demodified',
             'peptide_id',
             'protein_accession',
             'aa_start',
             'aa_end',
             'library'
            ],
        'representatives':
            [
             'library_id',
             'charge',
             'peaks',
             'psm_id',
             'variant_id',
             'peptide_id',
             'version_added',
             'score'
             ],
        'psm_provenance':
            [
             'file_descriptor',
             'nativeid',
             'dataset_id',
             'charge',
             'original_search',
             'workflow_name',
             'url',
             'variant_sequence'
            ],
        'provenance_representative':
            [
             'psm_id',
             'representative_id'
            ],
        'protein_datasets':
            [
             'dataset_id',
             'protein_id',
             'unique_peptides',
             'exon_unique',
             'splice_junction',
             'exon_mapped',
             'psms_shared',
             'psms_unique'
            ],
        'protein_statistics':
            [
             'protein',
             'library_id',
             'version',
             'unique_peptides',
             'exon_unique',
             'splice_junction',
             'exon_mapped',
             'peptides',
             'variants',
             'modifications',
             'datasets',
             'psms_shared',
             'psms_unique',
             'hupo_nonoverlapping'
            ]
    }[key]

def table_output(index,table,output_folder,output_rows, postprocess = ''):
    print("{}: {} rows".format(table,len(output_rows)));
    variant_output_file = os.path.join(output_folder,"{}_{}.tsv".format(index,table))
    variant_output_file_processed = variant_output_file
    if postprocess:
        variant_output_file_processed = os.path.join(output_folder,"{}{}_{}.tsv".format(index,postprocess,table))
    with open(variant_output_file, 'w') as w:
        for (i,row) in enumerate(output_rows):
            w.write("{}\n".format(",".join(["\"{}\"".format(r) for r in row])))


def mapping_to_peptides_and_mapping(input_mapping, library_id):
    output_peptides, output_mapping = [],[]
    for peptide, mapping_obj in input_mapping.items():
        output_peptides.append([
            peptide,
            mapping_obj['canonical_matches'],
            mapping_obj['gene_unique'],
            0,
            0,
            0,
            len(peptide)
        ])
        for protein,(start_aa, end_aa),exact in mapping_obj['mappings']:
            output_mapping.append([
                peptide,
                "#peptide_id:" + peptide,
                protein,
                start_aa,
                end_aa,
                library_id if exact else 0
            ])
    return output_peptides, output_mapping
