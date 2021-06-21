import argparse
import sys
import csv
from collections import defaultdict
import json
from pathlib import Path
from itertools import chain
import explorer_export
import read_mappings
from python_ms_utilities import mapping, resources

def arguments():
    parser = argparse.ArgumentParser(description='mzTab to list of peptides')
    parser.add_argument('-k','--comparison_pep', type = Path, help='Peptides to Compare')
    parser.add_argument('-f','--fasta', type = Path, help='Input FASTA')
    parser.add_argument('-a','--input_psms', type = Path, help='Input PSMs')
    parser.add_argument('-y','--input_psms_external', type = Path, help='Input PSMs (External)')
    parser.add_argument('-o','--output_psms', type = Path, help='Output PSMs')
    parser.add_argument('-p','--output_peptides', type = Path, help='Output Peptides')
    parser.add_argument('-r','--output_proteins', type = Path, help='Output Proteins')
    parser.add_argument('-b','--output_exons', type = Path, help='Output Exons')
    parser.add_argument('-c','--protein_coverage', type = Path, help='Added Protein Coverage')
    parser.add_argument('-d','--protein_coverage_external', type = Path, help='Added Protein Coverage (External)')
    parser.add_argument('-t','--cosine_cutoff', type = float, help='Cosine Cutoff')
    parser.add_argument('-l','--explained_intensity_cutoff', type = float, help='Explained Intensity Cutoff')
    parser.add_argument('-w','--nextprot_releases', type = Path, help='NextProt Releases')
    parser.add_argument('-m','--msv_to_pxd_mapping', type = Path, help='MSV to PXD Mapping')
    parser.add_argument('-g','--external_provenance', type = Path, help='Provenance from Library Creation')
    parser.add_argument('-v','--library_version', type = int, help='Library Version')
    parser.add_argument('-n','--library_name', type = int, help='Library Name')
    parser.add_argument('-x','--export_explorers', type = int, help='Export Explorer Tables (0/1)')
    parser.add_argument('-e','--explorers_output', type = Path, help='Tables for Explorers')

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

protein_type = lambda protein, proteome: 'TrEMBL' if proteome.proteins[protein].db == 'tr' else ('Canonical' if proteome.proteins[protein].iso == None else 'Isoform')

def msv_to_pxd(msv, msv_mapping):
    output_mapping = msv_mapping.get(msv,{}).get('px_accession')
    if not output_mapping:
        output_mapping = msv
    return output_mapping

def correct_usi(usi_input, msv_mapping):
    if '[' in usi_input or 'PXD' in usi_input:
        return usi_input
    else:
        split_usi = usi_input.split(':')
        split_usi[1] = msv_to_pxd(split_usi[1], msv_mapping)
        split_usi[5] = '/'.join([add_brackets(split_usi[5].split('/')[0]),split_usi[5].split('/')[1]])
        return ':'.join(split_usi)


def find_overlap(existing_peptides, new_peptides, protein_length, protein_pe, name, max_seq_output = 10):
    comparison_pos = set([sorted(positions, key = lambda x: x[0])[0] for positions in existing_peptides.values()])
    added_pos = set([sorted(positions, key = lambda x: x[0])[0] for positions in new_peptides.values()])

    comparison_peptides = set(["{} ({}-{})".format(peptide.replace('L','I'),sorted(positions, key = lambda x: x[0])[0][0],sorted(positions, key = lambda x: x[0])[0][1]) for peptide, positions in existing_peptides.items()])
    added_peptides = set(["{} ({}-{})".format(peptide.replace('L','I'),sorted(positions, key = lambda x: x[0])[0][0],sorted(positions, key = lambda x: x[0])[0][1]) for peptide, positions in new_peptides.items()])

    novel_peptides = added_peptides.difference(comparison_peptides)
    supporting_peptides = added_peptides.intersection(comparison_peptides)

    nonoverlapping_comparison_peptides = mapping.count_non_nested(comparison_pos.difference(added_pos))
    nonoverlapping_added_peptides = mapping.count_non_nested(added_pos.difference(comparison_pos))
    nonoverlapping_intersection = mapping.count_non_nested(comparison_pos.union(added_pos))
    nonoverlapping_all_previous = mapping.count_non_nested(comparison_pos)

    added_coverage = mapping.find_coverage(list(existing_peptides.values()), list(new_peptides.values()), protein_length)
    total_coverage = mapping.find_coverage([],list(existing_peptides.values()) + list(new_peptides.values()), protein_length)

    return {
        'supporting_peptides'+name:','.join(list(supporting_peptides)[:max_seq_output]) if len(supporting_peptides) >= 1 else ' ',
        'novel_peptides'+name:','.join(list(novel_peptides)[:max_seq_output]) if len(novel_peptides) >= 1 else ' ',
        'comparison_hpp'+name:nonoverlapping_all_previous,
        'new_hpp'+name:nonoverlapping_added_peptides,
        'combined_hpp'+name:nonoverlapping_intersection,
        # 'added_coverage'+name:added_coverage,
        'total_coverage'+name:total_coverage,
        'coverage_increase'+name:total_coverage/(total_coverage-added_coverage) if total_coverage-added_coverage != 0 else 1000,
        'promoted'+name:'Yes' if (nonoverlapping_all_previous < 2 and nonoverlapping_intersection >= 2) else 'No'
    }, [s.split(' ')[0] for s in supporting_peptides]

def protein_info(il_peptide, peptide_to_protein, all_proteins, added_proteins, proteome, comparison_seq, nextprot_pe):
    outdict = {}

    proteins = peptide_to_protein.get(il_peptide)

    if proteins:
        output_proteins = [protein for protein in proteins]
        cannonical_proteins = [protein for protein in proteins if proteome.proteins[protein].db == 'sp' and not proteome.proteins[protein].iso]
        noncannonical_proteins = [protein for protein in proteins if not (proteome.proteins[protein].db == 'sp' and not proteome.proteins[protein].iso)]
        output_genes = set([proteome.proteins[protein].gene if proteome.proteins[protein].gene else 'N/A' for protein in proteins])
        output_types = set([protein_type(protein, proteome) for protein in proteins])
    else:
        cannonical_proteins = []
        output_proteins = []
        output_genes = set()
        output_types = set()

    if proteins:
        if len(cannonical_proteins) > 0 and len(noncannonical_proteins) > 0:
            protein_str = ';'.join(cannonical_proteins) + " ###" + ';'.join(noncannonical_proteins)
        else:
            protein_str = ';'.join(output_proteins)
        outdict.update({
            'protein': protein_str,
            'gene': ';'.join(output_genes),
            'protein_type':';'.join(sorted(list(output_types))),
            'all_proteins': ';'.join(proteins)
        })
        if len(cannonical_proteins) == 1:
            protein_no_decoy = cannonical_proteins[0].replace('XXX_','')
            outdict.update({
                'pe': nextprot_pe.get(protein_no_decoy,''),
            })
        else:
            outdict.update({
                'pe': 'Multiple',
            })
    else:
        outdict.update({
            'protein': '',
            'pe': '',
        })
    if len(cannonical_proteins) == 1 and il_peptide in all_proteins.get(cannonical_proteins[0],{}):
        outdict['aa_start'],outdict['aa_end'] = all_proteins[cannonical_proteins[0]][il_peptide][0][0]
    elif len(output_proteins) == 1 and il_peptide in all_proteins.get(output_proteins[0],{}):
        outdict['aa_start'],outdict['aa_end'] = all_proteins[output_proteins[0]][il_peptide][0][0]
    else:
        outdict['aa_start'],outdict['aa_end'] = "N/A","N/A"

    if len(output_proteins) > 0 and len([g for g in output_genes if g != 'N/A']) <= 1 and len(cannonical_proteins) <= 1 and il_peptide in added_proteins[output_proteins[0]]:
        if il_peptide in comparison_seq:
            outdict['type'] = 'Matches existing evidence'
        else:
            outdict['type'] = 'New protein evidence'
        outdict['hpp_match'] = 'True'
    else:
        outdict['type'] = 'N/A'
        outdict['hpp_match'] = 'False'
    return outdict

def main():
    args = arguments()

    representative_per_precursor = {}

    proteome = mapping.add_decoys(mapping.read_uniprot(args.fasta))

    comparison_proteins = defaultdict(lambda: defaultdict(list))
    added_proteins = defaultdict(lambda: defaultdict(list))
    all_proteins = defaultdict(lambda: defaultdict(list))
    pep_mapping_info = {}
    comparison_proteins_w_synthetic = defaultdict(lambda: defaultdict(list))
    added_proteins_w_synthetic = defaultdict(lambda: defaultdict(list))
    has_synthetic = set()

    added_proteins_matching_synthetic = defaultdict(lambda: defaultdict(list))
    added_proteins_matching_synthetic_cosine = defaultdict(lambda: defaultdict(list))
    added_proteins_explained_intensity = defaultdict(lambda: defaultdict(list))
    added_proteins_isoform_unique = defaultdict(lambda: defaultdict(list))

    with open(args.msv_to_pxd_mapping) as json_file:
        msv_mapping = json.load(json_file)

    comparison_seq = set()

    with open(args.comparison_pep) as f:
        r = csv.DictReader(f, delimiter='\t')
        for l in r:
            comparison_seq.add(l['demodified'].replace('I','L'))
            comparison_proteins[l['protein']][l['demodified'].replace('I','L')].append((int(l['aa_start']),int(l['aa_end'])))

    for protein in comparison_proteins:
        for seq in comparison_proteins[protein].keys():
            if seq in has_synthetic:
                comparison_proteins_w_synthetic[protein][seq] = comparison_proteins[protein][seq]

    peptide_to_protein = defaultdict(list)
    peptide_to_protein_added = {}

    peptide_to_exon_map = defaultdict(list)

    nextprot_releases_pe = {}

    for nextprot_release in args.nextprot_releases.glob('*'):
        release_str,pe_dict = resources.get_nextprot_pe(nextprot_release)
        nextprot_releases_pe[release_str] = pe_dict

    latest_nextprot_release = sorted(list(nextprot_releases_pe.keys()))[-1]
    nextprot_pe = nextprot_releases_pe[latest_nextprot_release]


    for protein_coverage_file in chain(args.protein_coverage.glob('*'),args.protein_coverage_external.glob('*')):
        if(protein_coverage_file.is_file()):
            all_proteins,added_proteins,pep_mapping_info,peptide_to_exon_map = read_mappings.read_protein_coverage(protein_coverage_file,all_proteins,added_proteins,pep_mapping_info,peptide_to_exon_map,proteome)

    for protein in all_proteins:
        for peptide in all_proteins[protein].keys():
            peptide_to_protein[peptide.replace('I','L')].append(protein)

    for protein in added_proteins:
        if len(nextprot_pe) == 0 or protein in nextprot_pe:
            for peptide in added_proteins[protein].keys():
                peptide_to_protein_added[peptide.replace('I','L')] = protein
                if peptide.replace('I','L') in has_synthetic:
                    added_proteins_w_synthetic[protein][peptide] = added_proteins[protein][peptide]

    with open(args.output_psms,'w') as w:
        header = ['protein','protein_type','gene','all_proteins','decoy','pe','ms_evidence','filename','scan','sequence','sequence_unmodified','sequence_unmodified_il','charge','usi','score','modifications','pass','type','parent_mass','frag_tol','synthetic_filename','synthetic_scan','synthetic_usi','cosine','synthetic_match','explained_intensity','hpp_match','gene_unique','canonical_matches','all_proteins_w_coords','aa_start','aa_end', 'total_unique_exons_covered', 'exons_covered_no_junction', 'exon_junctions_covered', 'all_mapped_exons']

        o = csv.DictWriter(w, delimiter='\t',fieldnames = header, restval='N/A')
        o.writeheader()
        for input_psm in chain(args.input_psms.glob('*'),args.input_psms_external.glob('*')):
            with open(input_psm) as f:
                r = csv.DictReader(f, delimiter='\t')
                for l in r:
                    peptide = ''.join([p for p in l['sequence'] if p.isalpha()])
                    il_peptide = ''.join([p.replace('I','L') for p in l['sequence'] if p.isalpha()])

                    l['usi'] = correct_usi(l['usi'], msv_mapping)
                    l['synthetic_usi'] = correct_usi(l['synthetic_usi'], msv_mapping) if (l['synthetic_usi'] != 'N/A' and l['synthetic_usi'] != '') else 'N/A'
                    # l['sequence'] = l['sequence'] if '[' in l['sequence'] else add_brackets(l['sequence'])
                    l['sequence_unmodified'] = peptide
                    l['sequence_unmodified_il'] = il_peptide
                    l.update(pep_mapping_info.get(peptide,{}))
                    l.update(protein_info(il_peptide, peptide_to_protein, all_proteins, added_proteins, proteome, comparison_seq, nextprot_pe))
                    l.pop('mapped_proteins')
                    o.writerow(l)
                    sequence, charge = l['sequence'],l['charge']
                    if not (sequence, charge) in representative_per_precursor:
                        representative_per_precursor[(sequence, charge)] = l.copy()
                        precursor_representative = representative_per_precursor[(sequence, charge)]
                        precursor_representative['database_filename'] = l['filename']
                        precursor_representative['database_scan'] = l['scan']
                        precursor_representative['database_usi'] = l['usi']
                        precursor_representative['explained_intensity'] = float(l['explained_intensity'])
                        precursor_representative['cosine'] = float(l['cosine'])
                        precursor_representative['score'] = float(l['score'])
                        precursor_representative.pop('filename')
                        precursor_representative.pop('scan')
                        precursor_representative.pop('usi')

                    precursor_representative = representative_per_precursor[(sequence, charge)]

                    best_cosine = float(l['cosine']) >= precursor_representative['cosine']
                    best_score = float(l['score']) >= precursor_representative['score']
                    if best_score:
                        precursor_representative['database_filename'] = l['filename']
                        precursor_representative['database_scan'] = l['scan']
                        precursor_representative['database_usi'] = l['usi']
                        precursor_representative['score'] = float(l['score'])
                        precursor_representative['explained_intensity'] = float(l['explained_intensity'])
                    if best_cosine and float(l['cosine']) >= 0 and l['synthetic_usi'] != 'N/A':
                        precursor_representative['cosine_filename'] = l['filename']
                        precursor_representative['cosine_scan'] = l['scan']
                        precursor_representative['cosine_usi'] = l['usi']
                        precursor_representative['synthetic_filename'] = l['synthetic_filename']
                        precursor_representative['synthetic_scan'] = l['synthetic_scan']
                        precursor_representative['synthetic_usi'] = l['synthetic_usi']
                        precursor_representative['cosine'] = float(l['cosine'])
                        precursor_representative['explained_intensity'] = float(l['explained_intensity'])
                    if (best_score and best_cosine):
                        precursor_representative['cosine_score_match'] = 'Yes'
                    else:
                        precursor_representative['cosine_score_match'] = 'No'


    with open(args.output_peptides,'w') as w:
        header = ['protein','protein_type','gene','decoy','all_proteins','pe','ms_evidence','aa_total','database_filename','database_scan','database_usi','sequence','sequence_unmodified','sequence_unmodified_il','charge','score','modifications','pass','type','parent_mass','cosine_filename','cosine_scan','cosine_usi','synthetic_filename','synthetic_scan','synthetic_usi','cosine','synthetic_match','cosine_score_match','explained_intensity','hpp_match','gene_unique','canonical_matches','all_proteins_w_coords','aa_start','aa_end','frag_tol', 'total_unique_exons_covered', 'exons_covered_no_junction', 'exon_junctions_covered', 'all_mapped_exons']
        r = csv.DictWriter(w, delimiter = '\t', fieldnames = header, restval='N/A')
        r.writeheader()
        for (sequence, charge), best_psm in representative_per_precursor.items():
            sequence_nomod = ''.join([p for p in sequence if p.isalpha()])
            sequence_il = ''.join([p.replace('I','L') for p in sequence if p.isalpha()])
            best_psm.update(pep_mapping_info.get(sequence_nomod,{}))
            best_psm.update(protein_info(sequence_il, peptide_to_protein, all_proteins, added_proteins, proteome, comparison_seq, nextprot_pe))
            best_psm.pop('mapped_proteins')
            r.writerow(best_psm)

            proteins = peptide_to_protein.get(sequence_il,[])
            cannonical_proteins = [protein for protein in proteins if proteome.proteins[protein].db == 'sp' and not proteome.proteins[protein].iso]
            output_genes = set([proteome.proteins[protein].gene if proteome.proteins[protein].gene else 'N/A' for protein in proteins])

            if len([g for g in output_genes if g != 'N/A']) <= 1 and len(cannonical_proteins) <= 1:
                for protein in cannonical_proteins:
                    if sequence_il in added_proteins[protein] and best_psm['explained_intensity'] >= args.explained_intensity_cutoff:
                        if best_psm['cosine'] >= 0:
                            added_proteins_matching_synthetic[protein][sequence_il] = added_proteins[protein][sequence_il]
                            if best_psm['cosine'] >= args.cosine_cutoff:
                                added_proteins_matching_synthetic_cosine[protein][sequence_il] = added_proteins[protein][sequence_il]
                        added_proteins_explained_intensity[protein][sequence_il] = added_proteins[protein][sequence_il]
            if len(proteins) == 1:
                for protein in proteins:
                    if sequence_il in added_proteins[protein] and best_psm['explained_intensity'] >= args.explained_intensity_cutoff:
                        added_proteins_isoform_unique[protein][sequence_il] = added_proteins[protein][sequence_il]

    with open(args.output_proteins, 'w') as fo:

        overlap_keys = find_overlap({},{},0,0,"")[0].keys()

        fieldnames = [
            'protein',
            'protein_type',
            'gene',
            'pe',
            'aa_total',
            'ms_evidence',
            'fdr',
            'cosine_cutoff',
            'explained_intensity_cutoff'
        ] + ['_dyn_#neXtProt Release {}'.format(release) for release in sorted(nextprot_releases_pe.keys())]

        for overlap_suffix in ['','_w_synthetic','_w_synthetic_cosine','_just_current','_just_current_iso_unique']:
            fieldnames.extend(['{}{}'.format(k,overlap_suffix) for k in overlap_keys])

        w = csv.DictWriter(fo, delimiter = '\t', fieldnames = fieldnames)
        w.writeheader()

        for i,(protein,protein_entry) in enumerate(proteome.proteins.items()):

            is_canonical = protein_entry.db == 'sp' and not protein_entry.iso

            protein_dict = {
                'protein': protein,
                'protein_type':protein_type(protein, proteome),
                'pe': nextprot_pe.get(protein_entry.id,0) if is_canonical else 0,
                'aa_total':protein_entry.length,
                'gene':protein_entry.gene,
            }

            if is_canonical:
                protein_dict.update(find_overlap(comparison_proteins.get(protein,{}),added_proteins_explained_intensity.get(protein,{}),int(protein_dict['aa_total']),int(protein_dict['pe']),'')[0])
            else:
                protein_dict.update(find_overlap(comparison_proteins.get(protein,{}),added_proteins_isoform_unique.get(protein,{}),int(protein_dict['aa_total']),int(protein_dict['pe']),'')[0])

            protein_dict.update(find_overlap(comparison_proteins.get(protein,{}),added_proteins_matching_synthetic.get(protein,{}),int(protein_dict['aa_total']),int(protein_dict['pe']),'_w_synthetic')[0])
            protein_dict.update(find_overlap(comparison_proteins.get(protein,{}),added_proteins_matching_synthetic_cosine.get(protein,{}),int(protein_dict['aa_total']),int(protein_dict['pe']),'_w_synthetic_cosine')[0])
            protein_dict.update(find_overlap({},added_proteins_explained_intensity.get(protein,{}),int(protein_dict['aa_total']),int(protein_dict['pe']),'_just_current')[0])
            protein_dict.update(find_overlap({},added_proteins_isoform_unique.get(protein,{}),int(protein_dict['aa_total']),int(protein_dict['pe']),'_just_current_iso_unique')[0])

            protein_dict.update({
                'cosine_cutoff':args.cosine_cutoff,
                'explained_intensity_cutoff':args.explained_intensity_cutoff
            })

            for release, pe_dict in nextprot_releases_pe.items():
                print(release)
                protein_dict['_dyn_#neXtProt Release {}'.format(release)] = pe_dict.get(protein_entry.id,0) if is_canonical else 0
            w.writerow(protein_dict)

    with open(args.output_exons, 'w') as fo:

        fieldnames = [
            'chromosome',
            'start_location',
            'end_location',
            'gene',
            '#transcripts',
            '#peptides',
            '#peptides_novel',
            'coverage',
            'coverage_novel',
            '#junction_peptides',
            '#junction_peptides_novel'
        ]

        w = csv.DictWriter(fo, delimiter = '\t', fieldnames = fieldnames)
        w.writeheader()

        for exon, peptide_mappings in peptide_to_exon_map.items():
            chromosome, location, gene, transcripts = exon
            peptides = set([p[0] for p in peptide_mappings])
            junction_peptides = set([p[0] for p in peptide_mappings if len(p[0]) > 1])

            mapped_coordinates = []
            for peptide_mapping in peptide_mappings:
                for i,reference_mapping in enumerate(peptide_mapping[3]):
                    if reference_mapping == location:
                        mapped_coordinates.append(peptide_mapping[2][i])
            w.writerow({
                'chromosome':chromosome,
                'start_location':location[0],
                'end_location':location[1],
                'gene':gene,
                '#transcripts':len(transcripts),
                '#peptides':len(peptides),
                '#junction_peptides':len(junction_peptides),
                'coverage':mapping.find_coverage([mapped_coordinates],[],location[1],location[0]),
                '#peptides_novel':0,
                'coverage_novel':0,
                '#junction_peptides_novel':0
            })

    if args.export_explorers and args.export_explorers == 1:
        hupo_mapping = {}
        unique_mapping = {}
        for protein, peptides in added_proteins.items():
            if 'XXX_' not in protein:
                protein_dict = {
                    'protein': protein,
                    'pe': nextprot_pe.get(protein,0),
                    'aa_total':proteome.proteins.get(protein).length
                }
                hupo_mapping[protein] = find_overlap({},peptides,int(protein_dict['aa_total']),int(protein_dict['pe']),'')[0]['new_hpp']
        for protein, peptides in all_proteins.items():
            if 'XXX_' not in protein:
                unique_mapping[protein] = len(set([peptide for peptide, mappings in peptides.items() if mappings[0][1]]))
        explorer_export.output_for_explorer(args.explorers_output, pep_mapping_info, representative_per_precursor, args.external_provenance, hupo_mapping, unique_mapping, args.library_name, args.library_version)

if __name__ == '__main__':
    main()
