from pathlib import Path
from csv import DictReader, DictWriter
from collections import defaultdict

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
    variant_output_file = output_folder.joinpath("{}_{}.tsv".format(index,table))
    variant_output_file_processed = variant_output_file
    if postprocess:
        variant_output_file_processed = output_folder.joinpath("{}{}_{}.tsv".format(index,postprocess,table))
    with open(variant_output_file, 'w') as w:
        for (i,row) in enumerate(output_rows):
            w.write("{}\n".format(",".join(["\"{}\"".format(r) for r in row])))

def table_from_dataset_proteins(dataset_protein_dict, library_id, library_version):
    dataset_protein_rows = []
    for ((dataset,protein),d) in dataset_protein_dict.items():
        dataset_protein_rows.append([
            "#library_dataset_id:{}.{}".format(library_id,library_version),
            dataset,
            "#protein_id:" + protein,
            str(len(d['unique_peptides'])),
            str(len(d['exon_unique'].intersection(d['unique_peptides']))),
            str(len(d['splice_junction'].intersection(d['unique_peptides']))),
            str(len(d['exon_mapped'])),
            str(d['psms_shared']),
            str(d['psms_unique'])
        ])
    return dataset_protein_rows

def table_from_proteins(protein_dict, hupo_mapping, unique_mapping, library_id, library_version):
    protein_rows = []
    for (protein,d) in protein_dict.items():
        protein_rows.append([
            "#library_dataset_id:{}.{}".format(library_id,library_version),
            "#protein_id:" + protein,
            str(unique_mapping.get(protein,0)),
            str(len(d['exon_unique'].intersection(d['unique_peptides']))),
            str(len(d['splice_junction'].intersection(d['unique_peptides']))),
            str(len(d['exon_mapped'])),
            str(d['psms_shared']),
            str(d['psms_unique']),
            str(hupo_mapping.get(protein,0))
        ])
    return protein_rows

def mapping_to_peptides_and_mapping(input_mapping, library_id):
    peptides, mappings = [],[]
    for peptide, mapping_obj in input_mapping.items():
        peptides.append([
            peptide,
            mapping_obj['canonical_matches'],
            1 if mapping_obj['gene_unique'] else 2,
            1 if (mapping_obj['exon_junctions_covered'] + mapping_obj['exons_covered_no_junction']) == 1 else 0,
            mapping_obj['exon_junctions_covered'],
            mapping_obj['total_unique_exons_covered'],
            len(peptide)
        ])
        for protein,(start_aa,end_aa),il_ambiguous in mapping_obj['mapped_proteins']:
            mappings.append([
                peptide,
                "#peptide_id:" + peptide,
                "#protein_id:" + protein,
                start_aa,
                end_aa
            ])
    return peptides, mappings

def representatives_to_representatives_and_variants(representatives_table, library_id, library_version, task_for_filescan):
    representatives, variants = [], []
    for (sequence, charge), representative in representatives_table.items():
        spectrum_file_descriptor = "f.{}".format(representative["database_filename"].replace("jswertz/MSV000086369_hct116_symlinks", "MSV000086369/ccms_peak/RAW"))
        nativeid = "scan={}".format(representative["database_scan"])
        representatives.append([
            "#library_dataset_id:{}.{}".format(library_id,library_version),
            "#variant_id:{}.{}".format(sequence,charge),
            "#psm_id:{}.{}.{}".format(spectrum_file_descriptor,nativeid,task_for_filescan[(spectrum_file_descriptor,nativeid)]),
            representative['score']
        ])
        variants.append([
            representative['sequence'],
            representative['charge'],
            float(representative['parent_mass'])
        ])
    return representatives, variants

def update_provenance(input_provenance, input_representatives, input_mappings, library_id, library_version, output_provenance, output_provenance_representatives):
    proteins_stats = defaultdict(
        lambda : {
            'exon_unique':set(),
            'splice_junction':set(),
            'exon_mapped':set(),
            'unique_peptides':set(),
            'peptides':set(),
            'variants':set(),
            'modifications':set(),
            'datasets':set(),
            'psms_unique':0,
            'psms_shared':0,
            'hupo_nonoverlapping':0
        }
    )
    dataset_proteins = defaultdict(lambda : {'exon_unique':set(),'splice_junction':set(),'exon_mapped':set(), 'unique_peptides':set(), 'psms_unique':0, 'psms_shared':0})
    task_for_filescan = {}

    provenance_file = None
    provenance_generator = None

    if Path(input_provenance).is_file():
        print("External Provenance File")
        provenance_file = open(input_provenance)
        provenance_lines = DictReader(provenance_file, delimiter = '\t')
        get_sequence = lambda l: l['annotation']
        get_charge = lambda l: l['charge']
        # replace is because a symlink was used in construction of KB 2.0.1
        get_filename = lambda l: l['filename'].replace("jswertz/MSV000086369_hct116_symlinks", "MSV000086369/ccms_peak/RAW")
        get_scan = lambda l: l['scan']
        get_proteosafe_task = lambda l: l['proteosafe_task']
        get_workflow = lambda l, a: l.get('workflow',a)
        get_search_url = lambda l: l.get('search_url','https://proteomics2.ucsd.edu/ProteoSAFe/status.jsp?task={}'.format(get_proteosafe_task(l)))
    else:
        provenance_lines = input_provenance.values()
        get_sequence = lambda l: l['sequence']
        get_charge = lambda l: l['charge']
        get_filename = lambda l: l['filename']
        get_scan = lambda l: l['scan']
        get_proteosafe_task = lambda l: l['proteosafe_task']
        get_workflow = lambda l, a: l.get('workflow',a)
        get_search_url = lambda l: l.get('search_url','https://proteomics2.ucsd.edu/ProteoSAFe/status.jsp?task={}'.format(get_proteosafe_task(l)))

    with open(output_provenance_representatives, 'w') as w_provenance_representatives:
        for l in provenance_lines:
            if (get_sequence(l),get_charge(l)) in input_representatives:
                precursor_il = "".join([p.replace('I','L') for p in get_sequence(l)])
                sequence_il = "".join([p.replace('I','L') for p in get_sequence(l) if p.isalpha()])
                sequence_nomod = "".join([p for p in get_sequence(l) if p.isalpha()])
                dataset = get_filename(l).split('/')[0]
                mappings = input_mappings[sequence_nomod]
                # for peptide in peptide_il_map[il_peptide]:
                just_accessions = set([p[0] for p in mappings['mapped_proteins']])
                exon_unique = 1 if (mappings['exon_junctions_covered'] + mappings['exons_covered_no_junction']) == 1 else 0
                splice_junction = mappings['exon_junctions_covered']
                exon_mapped = mappings['total_unique_exons_covered']
                if len(just_accessions) > 0:
                    for protein in just_accessions:
                        proteins_stats[protein]['psms_shared'] += 1
                        dataset_proteins[(dataset,protein)]['psms_shared'] += 1
                        proteins_stats[protein]['peptides'].add(sequence_il)
                        proteins_stats[protein]['variants'].add(precursor_il)
                        proteins_stats[protein]['datasets'].add(dataset)
                        if len(just_accessions) == 1:
                            proteins_stats[protein]['psms_unique'] += 1
                            dataset_proteins[(dataset,protein)]['psms_unique'] += 1
                            dataset_proteins[(dataset,protein)]['unique_peptides'].add(sequence_il)
                            proteins_stats[protein]['unique_peptides'].add(sequence_il)
                        if int(exon_mapped) > 0:
                            dataset_proteins[(dataset,protein)]['exon_mapped'].add(sequence_il)
                            proteins_stats[protein]['exon_mapped'].add(sequence_il)
                        if int(splice_junction) > 0:
                            dataset_proteins[(dataset,protein)]['splice_junction'].add(sequence_il)
                            proteins_stats[protein]['splice_junction'].add(sequence_il)
                        if int(exon_unique) > 0:
                            dataset_proteins[(dataset,protein)]['exon_unique'].add(sequence_il)
                            proteins_stats[protein]['exon_unique'].add(sequence_il)
                # Jeremy 3/17/22: Commenting out the below file writes since we no longer need to
                # load this file into the merged ProteinExplorer/MassIVE search database schema
                #w_provenance.write(','.join(["\"{}\"".format(r) for r in [
                #    get_filename(l),
                #    'scan={}'.format(get_scan(l)),
                #    dataset,
                #    get_charge(l),
                #    get_proteosafe_task(l),
                #    get_workflow(l,'MSGF_PLUS').upper(),
                #    get_search_url(l),
                #    get_sequence(l)
                #]]) + '\n')
                w_provenance_representatives.write(','.join(["\"{}\"".format(r) for r in [
                    "#library_dataset_id:{}.{}".format(library_id,library_version),
                    "#variant_id:{}.{}".format(get_sequence(l),get_charge(l)),
                    '#psm_id:f.{}.scan={}.{}'.format(get_filename(l),get_scan(l),get_proteosafe_task(l))
                ]]) + '\n')
                representative = input_representatives[(get_sequence(l),get_charge(l))]
                if get_filename(l) == representative["database_filename"] and get_scan(l) == representative["database_scan"]:
                    task_for_filescan[("f.{}".format(get_filename(l)),"scan={}".format(get_scan(l)))] = get_proteosafe_task(l)
    if Path(input_provenance).is_file():
        provenance_file.close()
    return dataset_proteins, proteins_stats, task_for_filescan


def output_for_explorer(output_folder, input_mappings, input_representatives, input_provenance, hupo_mapping, unique_mapping, library_id, library_version):
    # need to stream provenance
    output_provenance = output_folder.joinpath("{}_{}.tsv".format(8,'psm_provenance'))
    output_provenance_representatives = output_folder.joinpath("{}_{}.tsv".format(5,'library_variant_psms'))
    dataset_protein_dict, proteins_stats_dict, task_for_filescan_dict = update_provenance(input_provenance, input_representatives, input_mappings, library_id, library_version, output_provenance, output_provenance_representatives)

    peptides, mappings = mapping_to_peptides_and_mapping(input_mappings, library_id)
    representatives, variants = representatives_to_representatives_and_variants(input_representatives, library_id, library_version, task_for_filescan_dict)

    dataset_proteins = table_from_dataset_proteins(dataset_protein_dict, library_id, library_version)
    protein_stats = table_from_proteins(proteins_stats_dict, hupo_mapping, unique_mapping, library_id, library_version)

    table_output(1, 'peptides', output_folder, peptides)
    table_output(2, 'variants', output_folder, variants)
    table_output(3, 'peptide_proteins', output_folder, mappings)
    table_output(4, 'library_variants', output_folder, representatives, 'a')
    # table_output(5, 'library_variant_psms', output_folder, provenance_representatives, 'b')
    table_output(6, 'dataset_proteins', output_folder, protein_stats)
    table_output(7, 'library_dataset_proteins', output_folder, dataset_proteins)
    # table_output(8, 'psm_provenance', args.output_folder, provenance)
