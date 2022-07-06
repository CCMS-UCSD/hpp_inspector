import argparse
import re
import sys
import csv
from collections import defaultdict, namedtuple
import json
from pathlib import Path
from itertools import chain
import explorer_export
import read_mappings
from python_ms_utilities import mapping, resources, fdr, proteosafe
import pandas as pd
from datetime import datetime
import requests

def arguments():
    parser = argparse.ArgumentParser(description='mzTab to list of peptides')
    parser.add_argument('--comparison_pep', type = Path, help='Peptides to Compare')
    parser.add_argument('--proteome_fasta', type = Path, help='Input Proteome FASTA')
    parser.add_argument('--contaminants_fasta', type = Path, help='Input Contaminants FASTA')
    parser.add_argument('--input_psms', type = Path, help='Input PSMs')
    parser.add_argument('--input_psms_external', type = Path, help='Input PSMs (External)')
    parser.add_argument('--input_peptides', type = Path, help='Input PSMs (External)')
    parser.add_argument('--output_psms_flag', type = str, help='Output PSMs Flag')
    parser.add_argument('--output_psms', type = Path, help='Output PSMs')
    parser.add_argument('--output_peptides', type = Path, help='Output Peptides')
    parser.add_argument('--output_proteins', type = Path, help='Output Proteins')
    parser.add_argument('--output_exons', type = Path, help='Output Exons')
    parser.add_argument('--output_mappings', type = Path, help='Output Mappings')
    parser.add_argument('--output_dataset_proteins_hpp', type = Path, help='Output Dataset Proteins (HPP)')
    parser.add_argument('--output_dataset_proteins_all', type = Path, help='Output Dataset Proteins (All)')
    parser.add_argument('--output_task_proteins_hpp', type = Path, help='Output Task Proteins (HPP)')
    parser.add_argument('--output_task_proteins_all', type = Path, help='Output Task Proteins (All)')
    parser.add_argument('--protein_coverage', type = Path, help='Added Protein Coverage')
    parser.add_argument('--protein_coverage_external', type = Path, help='Added Protein Coverage (External)')
    parser.add_argument('--cosine_cutoff', type = float, help='Cosine Cutoff')
    parser.add_argument('--explained_intensity_cutoff', type = float, help='Explained Intensity Cutoff')
    parser.add_argument('--annotated_ions_cutoff', type = float, help='Annotated Ion Cutoff')
    parser.add_argument('--precursor_fdr', type = float, help='Precursor FDR')
    parser.add_argument('--picked_protein_fdr', type = float, help='Canonical Protein FDR')
    parser.add_argument('--hpp_protein_fdr', type = float, help='HPP Protein FDR')
    parser.add_argument('--picked_protein_fdr_comparison', type = float, help='Canonical Protein FDR')
    parser.add_argument('--hpp_protein_fdr_comparison', type = float, help='HPP Protein FDR')
    parser.add_argument('--filter_rows_fdr', type = int, help='Filter Rows By FDR')
    parser.add_argument('--main_fdr', type = str, help='Main FDR Type')
    parser.add_argument('--leftover_fdr', type = str, help='Leftover FDR Type')
    parser.add_argument('--nextprot_releases', type = Path, help='NextProt Releases')
    parser.add_argument('--msv_to_pxd_mapping', type = Path, help='MSV to PXD Mapping')
    parser.add_argument('--external_provenance', type = Path, help='Provenance from Library Creation')
    parser.add_argument('--library_version', type = int, help='Library Version')
    parser.add_argument('--library_name', type = int, help='Library Name')
    parser.add_argument('--export_explorers', type = int, help='Export Explorer Tables (0/1)')
    parser.add_argument('--explorers_output', type = Path, help='Tables for Explorers')
    if len(sys.argv) < 4:
        parser.print_help()
        sys.exit(1)
    return parser.parse_args()

aa_weights =    [71,  156, 114, 115, 103, 129, 128, 57,  137, 113, 113, 128, 131, 147, 97,  87,  101, 186, 163, 99]
aa_characters = ['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

aa_dict = dict(zip(aa_characters,aa_weights))

def unit_delta(aa1, aa2):
    return (aa_dict.get(aa2,0) - aa_dict.get(aa1,0))

theoretical_mass = lambda xs: sum([aa_dict[x] for x in xs]) + 18.010564686

no_mod_il = lambda pep: ''.join([p.replace('I','L') for p in pep if p.isalpha()])

def theoretical_mz(sequence,charge):
    aa = ''.join([a for a in sequence if a.isalpha()])
    mods = ''.join([m for m in sequence if not m.isalpha()])
    if len(mods) > 0:
        mods = eval(mods)
    else:
        mods = 0
    return (theoretical_mass(aa) + mods + (int(charge)*1.007276035))/int(charge)


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

protein_type = lambda protein, proteome: 'Contaminant' if proteome.proteins[protein].db == 'con' else ('TrEMBL' if proteome.proteins[protein].db == 'tr' else ('Canonical' if proteome.proteins[protein].iso == None else 'Isoform'))

def msv_to_pxd(msv, msv_mapping):
    output_mapping = msv_mapping.get(msv,{}).get('px_accession')
    if not output_mapping:
        output_mapping = msv
    return output_mapping

def correct_usi(usi_input, msv_mapping):
    try:
        if '[' in usi_input or 'PXD' in usi_input:
            return usi_input
        else:
            split_usi = usi_input.split(':')
            split_usi[1] = msv_to_pxd(split_usi[1], msv_mapping)
            split_usi[5] = '/'.join([add_brackets(split_usi[5].split('/')[0]),split_usi[5].split('/')[1]])
            return ':'.join(split_usi)
    except:
        return ''

def find_overlap(existing_peptides, new_peptides, protein_length, protein_pe, name, max_seq_output = 10, pass_this_picked_fdr = False, pass_this_hpp_fdr = False, pass_comparison_picked_fdr = False, pass_comparison_hpp_fdr = False):
    comparison_pos = set([sorted(positions, key = lambda x: x[0])[0] for positions in existing_peptides.values()])
    added_pos = set([sorted(positions, key = lambda x: x[0])[0] for positions in new_peptides.values()])

    comparison_peptides = set(["{} ({}-{})".format(peptide.replace('L','I'),sorted(positions, key = lambda x: x[0])[0][0],sorted(positions, key = lambda x: x[0])[0][1]) for peptide, positions in existing_peptides.items()])
    added_peptides = set(["{} ({}-{})".format(peptide.replace('L','I'),sorted(positions, key = lambda x: x[0])[0][0],sorted(positions, key = lambda x: x[0])[0][1]) for peptide, positions in new_peptides.items()])

    novel_peptides = added_peptides.difference(comparison_peptides)
    supporting_peptides = added_peptides.intersection(comparison_peptides)

    nonoverlapping_union = mapping.count_non_nested(added_pos.union(comparison_pos))
    nonoverlapping_added_peptides = mapping.count_non_nested(added_pos.difference(comparison_pos))
    nonoverlapping_all_previous = mapping.count_non_nested(comparison_pos)
    nonoverlapping_all_current = mapping.count_non_nested(added_pos)

    nonoverlapping_fdr_controlled = 0

    if (pass_comparison_hpp_fdr or pass_this_hpp_fdr):
        nonoverlapping_fdr_controlled = nonoverlapping_union
    else:
        nonoverlapping_fdr_controlled = 1 if nonoverlapping_union > 0 else 0


    nonoverlapping_all_previous_allfdr = 0
    nonoverlapping_all_previous_hppfdr = 0


    if pass_comparison_hpp_fdr:
        nonoverlapping_all_previous_hppfdr = nonoverlapping_all_previous
    elif pass_comparison_picked_fdr:
        nonoverlapping_all_previous_allfdr = nonoverlapping_all_previous

    added_coverage = mapping.find_coverage(list(existing_peptides.values()), list(new_peptides.values()), protein_length)
    total_coverage = mapping.find_coverage([],list(existing_peptides.values()) + list(new_peptides.values()), protein_length)

    return {
        'supporting_peptides'+name:','.join(list(supporting_peptides)[:max_seq_output]) if len(supporting_peptides) >= 1 else ' ',
        'novel_peptides'+name:','.join(list(novel_peptides)[:max_seq_output]) if len(novel_peptides) >= 1 else ' ',
        'comparison_hpp_hppfdr'+name:nonoverlapping_all_previous_hppfdr,
        'comparison_hpp_allfdr'+name:nonoverlapping_all_previous_allfdr,
        'new_hpp'+name:nonoverlapping_added_peptides,
        'just_added_hpp'+name:nonoverlapping_all_current,
        'combined_fdr_hpp'+name:nonoverlapping_fdr_controlled,
        'combined_hpp'+name:nonoverlapping_union,
        'total_coverage'+name:total_coverage,
        'coverage_increase'+name:total_coverage/(total_coverage-added_coverage) if total_coverage-added_coverage != 0 else 1000,
    }, [s.split(' ')[0] for s in supporting_peptides]

def main():
    args = arguments()

    representative_per_precursor = {}
    variant_to_best_precursor = {}

    proteome = mapping.add_decoys(mapping.merge_proteomes([mapping.read_uniprot(args.proteome_fasta),mapping.read_fasta(args.contaminants_fasta)]))

    all_datasets = set()
    datasets_per_sequence= defaultdict(set)
    all_tasks = set()
    tasks_per_sequence= defaultdict(set)

    SeqOccurances = namedtuple('SeqOccurances',['match','synthetic_match','synthetic_match_cosine'])
    SeqInfo = namedtuple('SeqInfo',['hpp','isoform_unique','added','comparison'])

    protein_mapping = defaultdict(lambda: defaultdict(set))
    sequences_found = defaultdict(lambda: SeqInfo(False,False,SeqOccurances(False,False,False),SeqOccurances(False,False,False)))

    pep_mapping_info = {}

    comparison_picked_fdr = {}
    comparison_hpp_fdr = {}

    frequency = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))

    with open(args.msv_to_pxd_mapping) as json_file:
        msv_mapping = json.load(json_file)
        pxd_mapping = {}
        for msv, msv_map in msv_mapping.items():
            pxd = msv_map.get('px_accession')
            if pxd:
                pxd_mapping[pxd] = msv


    with open(args.comparison_pep) as f:
        r = csv.DictReader(f, delimiter='\t')
        for l in r:
            il_peptide = l['demodified'].replace('I','L')
            has_synthetic = False
            has_synthetic_cosine = False
            if l['protein'] in proteome.proteins:
                
                if proteome.proteins[l['protein']].db == 'sp' and not proteome.proteins[l['protein']].iso:
                    #peptides that are not uniquely matching will have both FDRs as 1
                    comparison_picked_fdr[l['protein']] = min(float(l['all_protein_fdr']),comparison_picked_fdr.get(l['protein'],1))
                    comparison_hpp_fdr[l['protein']] = min(float(l['hpp_protein_fdr']),comparison_hpp_fdr.get(l['protein'],1))

                protein_mapping[l['protein']][il_peptide].add((int(l['aa_start']),int(l['aa_end']),None,None,None,None))
                comparison = sequences_found[il_peptide].comparison
                if comparison != SeqOccurances(True,True,True):
                    if float(l.get('synthetic_cosine',-1)) >= 0:
                        has_synthetic = True
                        if float(l['synthetic_cosine']) > args.cosine_cutoff:
                            has_synthetic_cosine = True
                    comparison = SeqOccurances(True, comparison.synthetic_match or has_synthetic, comparison.synthetic_match_cosine or has_synthetic_cosine)
                    sequences_found[il_peptide] = sequences_found[il_peptide]._replace(comparison = comparison)
                if l['is_hpp'] == 'True':
                    sequences_found[il_peptide] = sequences_found[il_peptide]._replace(hpp = True)

    peptide_to_protein = defaultdict(list)
    peptide_to_exon_map = defaultdict(list)

    nextprot_releases_pe = {}

    for nextprot_release in args.nextprot_releases.glob('*'):
        release_str,pe_dict = resources.get_nextprot_pe(nextprot_release)
        nextprot_releases_pe[release_str] = pe_dict

    print("Loaded PEs")

    latest_nextprot_release = sorted(list(nextprot_releases_pe.keys()))[-1]
    nextprot_pe = nextprot_releases_pe[latest_nextprot_release]

    def update_precursor_representative(l,from_psm = True):
        datasets = set([d for d in l.get('datasets','').split(';') if d != ''])
        tasks = set([t for t in l.get('tasks','').split(';') if t != ''])
        sequence, charge = l['sequence'],l['charge']
        precursor_theoretical_mz = theoretical_mz(sequence, charge)
        aa_seq = ''.join([a for a in sequence if a.isalpha()])

        variant = (aa_seq, charge, int(precursor_theoretical_mz))
        update_peptidoform = False

        if not variant in variant_to_best_precursor:
            variant_to_best_precursor[variant] = (sequence, charge)
            representative_per_precursor[(sequence, charge)] = l.copy()
            precursor_representative = representative_per_precursor[(sequence, charge)]
            precursor_representative['datasets'] = datasets
            precursor_representative['tasks'] = tasks
            precursor_representative['parent_mass'] = precursor_theoretical_mz
            if from_psm:
                precursor_representative['database_filename'] = l['filename']
                precursor_representative['database_scan'] = l['scan']
                precursor_representative['database_usi'] = l['usi']
                precursor_representative['explained_intensity'] = float(l['explained_intensity'])
                precursor_representative['cosine'] = float(l['cosine'])
                precursor_representative['score'] = float(l['score'])
                precursor_representative.pop('filename')
                precursor_representative.pop('scan')
                precursor_representative.pop('usi')

        precursor_representative = representative_per_precursor[variant_to_best_precursor[variant]]

        precursor_representative['datasets'] |= datasets
        precursor_representative['tasks'] |= tasks

        best_cosine = float(l['cosine']) >= float(precursor_representative['cosine'])
        best_score = float(l['score']) >= float(precursor_representative['score'])

        this_pass_ei = float(l['explained_intensity']) >= args.explained_intensity_cutoff
        this_pass_cos = float(l['cosine']) >= args.cosine_cutoff
        this_pass_by = int(l['matched_ions']) >= args.annotated_ions_cutoff
        best_pass_ei = float(precursor_representative['explained_intensity']) >= args.explained_intensity_cutoff
        best_pass_cos = float(precursor_representative['cosine']) >= args.cosine_cutoff
        best_pass_by = int(precursor_representative['matched_ions']) >= args.annotated_ions_cutoff

        # we want to find the spectrum with the best score,
        # but sometimes the highest scoring spectrum does not pass the filters

        potential_psm_loss = (
            ((not this_pass_ei and not this_pass_cos) or not this_pass_by) and
            ((best_pass_ei or best_pass_cos) and best_pass_by)
        )

        potential_psm_gain = (
            ((not best_pass_ei and not best_pass_cos) or not best_pass_by) and
            ((this_pass_ei or this_pass_cos) and this_pass_by)
        )

        if potential_psm_gain or (best_score and not potential_psm_loss):
            precursor_representative['database_filename'] = l['filename'] if from_psm else l['database_filename']
            precursor_representative['database_scan'] = l['scan'] if from_psm else l['database_scan']
            precursor_representative['database_usi'] = l['usi'] if from_psm else l['database_usi']
            precursor_representative['score'] = float(l['score'])
            update_peptidoform = True
            #consider best EI And matched ions over all representatives
            if potential_psm_gain or float(l['explained_intensity']) > float(precursor_representative.get('explained_intensity',0.0)):
                precursor_representative['explained_intensity'] = l['explained_intensity']
                precursor_representative['matched_ions'] = l['matched_ions']

        if best_cosine and float(l['cosine']) >= 0 and l['synthetic_usi'] != 'N/A':
            precursor_representative['cosine_filename'] = l['filename'] if from_psm else l['cosine_filename']
            precursor_representative['cosine_scan'] = l['scan'] if from_psm else l['cosine_scan']
            precursor_representative['cosine_usi'] = l['usi'] if from_psm else l['cosine_usi']
            precursor_representative['synthetic_filename'] = l['synthetic_filename']
            precursor_representative['synthetic_scan'] = l['synthetic_scan']
            precursor_representative['synthetic_usi'] = l['synthetic_usi']
            precursor_representative['cosine'] = float(l['cosine'])

        if (best_score and best_cosine):
            precursor_representative['cosine_score_match'] = 'Yes'
        else:
            precursor_representative['cosine_score_match'] = 'No'

        if update_peptidoform:
            representative_per_precursor[(sequence, charge)] = representative_per_precursor.pop(variant_to_best_precursor[variant])
            variant_to_best_precursor[variant] = (sequence, charge)


    def update_mappings(protein_coverage_file,update_precursor_representatives):
        if(protein_coverage_file.is_file()):
            print("{}: Loading {} ({} cumulative peptides @ {} peptides/second)".format(datetime.now().strftime("%H:%M:%S"),protein_coverage_file,len(pep_mapping_info),len(pep_mapping_info)/(1+(datetime.now()-start_time).seconds)))
            protein_mapping_added,pep_mapping_info_added,peptide_to_exon_map_added,output_peptides,_,_ = read_mappings.read_protein_coverage(protein_coverage_file,set(pep_mapping_info.keys()),proteome)
            pep_mapping_info.update(pep_mapping_info_added)
            peptide_to_exon_map.update(peptide_to_exon_map_added)
            for protein, peptide_mapping in protein_mapping_added.items():
                protein_mapping[protein].update(peptide_mapping)
            if update_precursor_representatives:
                for l in output_peptides:
                    update_precursor_representative(l,False)


    start_time = datetime.now()
    for protein_coverage_file in args.protein_coverage.glob('*'):
        update_mappings(protein_coverage_file,False)
            
    for protein_coverage_file in args.protein_coverage_external.glob('*'):
        update_mappings(protein_coverage_file,True)
            
    print("Loaded Protein Coverage")

    for protein in protein_mapping:
        for peptide in protein_mapping[protein].keys():
            peptide_to_protein[peptide.replace('I','L')].append(protein)

    for peptide, pep_mapping in pep_mapping_info.items():
        sequences_found[peptide.replace('I','L')] = sequences_found[peptide.replace('I','L')]._replace(hpp = pep_mapping['hpp'])
    
    print("About to load PSMs")


    def protein_info(peptide, peptide_to_protein, protein_mappings, sequences_found, proteome, nextprot_pe, spectrum_quality_pass):
        outdict = {}

        il_peptide = peptide.replace('I','L')
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
            if len(cannonical_proteins) > 0 and len(noncannonical_proteins) > 0:
                protein_str_full = ';'.join([mapping.protein_to_full_uniprot_string(proteome.proteins[p]) for p in cannonical_proteins]) + " ###" + ';'.join([mapping.protein_to_full_uniprot_string(proteome.proteins[p]) for p in noncannonical_proteins])
            else:
                protein_str_full = ';'.join([mapping.protein_to_full_uniprot_string(proteome.proteins[p]) for p in output_proteins])
            outdict.update({
                'protein': protein_str,
                'protein_full': protein_str_full,
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
        if len(cannonical_proteins) == 1 and il_peptide in protein_mappings.get(cannonical_proteins[0],{}):
            outdict['aa_start'],outdict['aa_end'],_,_,_,_ = next(iter(protein_mappings[cannonical_proteins[0]][il_peptide]))
        elif len(output_proteins) == 1 and il_peptide in protein_mappings.get(output_proteins[0],{}):
            outdict['aa_start'],outdict['aa_end'],_,_,_,_ = next(iter(protein_mappings[output_proteins[0]][il_peptide]))
        else:
            outdict['aa_start'],outdict['aa_end'] = "N/A","N/A"

        if sequences_found[il_peptide].hpp and len(output_proteins) > 0 and len(cannonical_proteins) <= 1 and il_peptide in protein_mappings[output_proteins[0]]:
            if sequences_found[il_peptide].comparison.match:
                if comparison_picked_fdr.get(output_proteins[0],1) <= 0.01:
                    outdict['type'] = 'Matches existing evidence'
                else:
                    outdict['type'] = 'New protein evidence (hint in reference)'
            else:
                outdict['type'] = 'New protein evidence'
            if spectrum_quality_pass:
                outdict['hpp_match'] = 'Yes'
            else:
                outdict['hpp_match'] = 'No - failed quality thresholds'
        else:
            outdict['type'] = 'Not HPP compliant'
            if spectrum_quality_pass:
                outdict['hpp_match'] = 'No - 2+ SAAV protein matches'
            else:
                outdict['hpp_match'] = 'No - failed quality thresholds'
        return outdict, output_proteins

    row_pass_filters = lambda best_psm: (float(best_psm['explained_intensity']) >= args.explained_intensity_cutoff or float(best_psm['cosine']) >= args.cosine_cutoff) and int(best_psm['matched_ions']) >= args.annotated_ions_cutoff

    if args.output_psms:
        all_psms_with_score = []
        all_psm_rows = []
        with open(args.output_psms,'w') as w:
            header = ['protein','protein_full','psm_fdr','protein_type','gene','all_proteins','decoy','pe','ms_evidence','filename','scan','sequence','sequence_unmodified','sequence_unmodified_il','charge','usi','score','modifications','pass','type','parent_mass','frag_tol','synthetic_filename','synthetic_scan','synthetic_usi','cosine','synthetic_match','explained_intensity','matched_ions','hpp_match','gene_unique','canonical_matches','all_proteins_w_coords','aa_start','aa_end', 'total_unique_exons_covered', 'exons_covered_no_junction', 'exon_junctions_covered', 'all_mapped_exons','datasets','tasks']
            o = csv.DictWriter(w, delimiter='\t',fieldnames = header, restval='N/A', extrasaction='ignore')
            o.writeheader()
            for input_psm in args.input_psms.glob('*'):
                with open(input_psm) as f:
                    r = csv.DictReader(f, delimiter='\t')
                    for l in r:
                        peptide = ''.join([p for p in l['sequence'] if p.isalpha()])
                        il_peptide = ''.join([p.replace('I','L') for p in l['sequence'] if p.isalpha()])

                        l['usi'] = correct_usi(l['usi'], msv_mapping)
                        l['datasets'] = l['usi'].split(':')[1] if l['usi'] != '' else 'User spectrum files'
                        # if float(l['explained_intensity']) >= args.explained_intensity_cutoff:
                        #     sequences_per_dataset[l['usi'].split(':')[1]].add(il_peptide)
                        l['synthetic_usi'] = correct_usi(l['synthetic_usi'], msv_mapping) if (l['synthetic_usi'] != 'N/A' and l['synthetic_usi'] != '') else 'N/A'
                        # l['sequence'] = l['sequence'] if '[' in l['sequence'] else add_brackets(l['sequence'])
                        l['sequence_unmodified'] = peptide
                        l['sequence_unmodified_il'] = il_peptide

                        l.update(pep_mapping_info.get(peptide,{}))

                        protein_info_dict, output_proteins = protein_info(peptide, peptide_to_protein, protein_mapping, sequences_found, proteome, nextprot_pe, row_pass_filters(l))
                        for protein in output_proteins:
                            aa_start,aa_end,_,_,_,_ = next(iter(protein_mapping[protein][il_peptide]))
                            frequency[protein][(aa_start,aa_end)][(l['sequence'],l['charge'])] += 1
                        l.update(protein_info_dict)
                        proteins = l['protein'].split(' ###')[0].split(';')
                        # PSM-level FDR was inefficient at this scale - need to rethink
                        #if row_pass_filters(l):
                        all_targets = [p for p in proteins if 'XXX_' not in p]
                        all_psms_with_score.append(fdr.ScoredElement(l['usi'],len(all_targets)==0,l['score']))
                        l.pop('mapped_proteins')
                        l.pop('hpp')
                        l.pop('len')
                        all_psm_rows.append(l)

            psm_fdr = fdr.calculate_fdr(all_psms_with_score)

            for l in all_psm_rows:
                l['psm_fdr'] = psm_fdr.get(l['usi'],1)
                if args.output_psms_flag == "1" or (args.output_psms_flag == "0.5" and row_pass_filters(l)):
                    o.writerow(l)
                if l['psm_fdr'] <= 0.01:
                    update_precursor_representative(l)


    print("About to calculate precursor FDR")

    all_hpp_precursors = []

    for (sequence, charge), best_psm in representative_per_precursor.items():
        all_targets = [p for ps in best_psm['protein'].split(' ###') for p in ps.split(';') if 'XXX_' not in p]
        if best_psm['protein'] != '':
            all_hpp_precursors.append(fdr.ScoredElement((sequence, charge),len(all_targets)==0,best_psm['score']))
    
    precursor_fdr = {}
    if len(all_hpp_precursors) > 0:
        precursor_fdr = fdr.calculate_fdr(all_hpp_precursors)

    precursors_per_protein_all = defaultdict(lambda: defaultdict(list))
    precursors_per_protein_hpp = defaultdict(lambda: defaultdict(list))
    precursors_per_protein_non_unique = defaultdict(lambda: defaultdict(list))


    def output_protein_level_results(best_psm):

        sequence = best_psm['sequence_unmodified']
        sequence_il = best_psm['sequence_unmodified_il']

        proteins = [p for p in best_psm['protein'].split(' ###')[0].split(';') if p != '']
        if row_pass_filters(best_psm):
            if float(best_psm['precursor_fdr']) <= args.precursor_fdr and len(proteins) == 1 and ('Canonical' in best_psm.get('protein_type','') or 'Contaminant' in best_psm.get('protein_type','')):
                pos = (int(best_psm['aa_start']),int(best_psm['aa_end']))
                if best_psm.get('hpp_match','') == 'Yes':
                    precursors_per_protein_hpp[proteins[0]][pos].append((float(best_psm['score']),theoretical_mz(best_psm['sequence'],int(best_psm['charge']))))
                precursors_per_protein_all[proteins[0]][pos].append((float(best_psm['score']),theoretical_mz(best_psm['sequence'],int(best_psm['charge']))))
            for protein in proteins:
                precursors_per_protein_non_unique[protein][sequence_il].append((float(best_psm['score']),theoretical_mz(best_psm['sequence'],int(best_psm['charge']))))

        proteins = peptide_to_protein.get(sequence_il,[])
        cannonical_proteins = [protein for protein in proteins if (proteome.proteins[protein].db == 'sp' and not proteome.proteins[protein].iso) or proteome.proteins[protein].db == 'con']
        output_genes = set([proteome.proteins[protein].gene if proteome.proteins[protein].gene else 'N/A' for protein in proteins])

        if len(cannonical_proteins) <= 1 and best_psm.get('hpp_match','') == 'Yes':
            is_hpp = pep_mapping_info[sequence]['hpp']
            match = False
            has_synthetic = False
            has_synthetic_cosine = False
            is_isoform_unique = False
            if len(cannonical_proteins) == 1:
                if row_pass_filters(best_psm):
                    if float(best_psm['cosine']) >= 0:
                        has_synthetic = True
                        if float(best_psm['cosine']) >= args.cosine_cutoff:
                            has_synthetic_cosine = True
                    match = True
                    for dataset in best_psm['datasets']:
                        all_datasets.add(dataset)
                        datasets_per_sequence[sequence_il].add(dataset)
                    for task in best_psm['tasks']:
                        all_tasks.add(task)
                        tasks_per_sequence[sequence_il].add(task)
            if len(proteins) == 1:
                if row_pass_filters(best_psm):
                    is_isoform_unique = True
            added = sequences_found[sequence_il].added
            sequences_found[sequence_il] = sequences_found[sequence_il]._replace(
                hpp=is_hpp,
                isoform_unique=is_isoform_unique,
                added = SeqOccurances(
                    added.match or match,
                    added.synthetic_match or has_synthetic,
                    added.synthetic_match_cosine or has_synthetic_cosine
                    )
                )
    print("About to output sequences")

    all_precursors = []

    if args.output_peptides:
        for (sequence, charge), best_psm in representative_per_precursor.items():
            sequence_nomod = ''.join([p for p in sequence if p.isalpha()])
            best_psm.update(pep_mapping_info.get(sequence_nomod,{}))
            best_psm.update(protein_info(sequence_nomod, peptide_to_protein, protein_mapping, sequences_found, proteome, nextprot_pe, row_pass_filters(best_psm))[0])
            best_psm.pop('mapped_proteins')
            best_psm.pop('hpp')
            best_psm.pop('len')
            best_psm['precursor_fdr'] = min(precursor_fdr.get((sequence, charge),1),1)
            best_psm['psm_fdr'] = -1
            best_psm['synthetic_sequence'] = sequence.replace('+229.163','').replace('+229.162932','')
            output_protein_level_results(best_psm)
            best_psm['datasets'] = ';'.join(best_psm['datasets'])
            best_psm['tasks'] = ';'.join(best_psm['tasks'])
            all_precursors.append(best_psm)

    elif args.input_peptides:
        with open(args.input_peptides) as f:
            r = csv.DictReader(f, delimiter='\t')
            for best_psm in r:
                output_protein_level_results(best_psm)

    hpp_protein_w_scores = []
    hpp_score_dict = {}

    leftover_protein_w_scores = []
    leftover_score_dict = {}

    picked_protein_w_scores = []
    picked_score_dict = {}

    common_protein_w_scores = []
    common_score_dict = {}

    seen_picked = set()

    def greedy_sequence_precursor_score(precursor_list, mz_distance = 2.5):
        used_precursors = []
        for precursor in sorted(precursor_list, key = lambda x: x[0], reverse = True):
            found = False
            for seen_precursor in used_precursors:
                if abs(precursor[1]-seen_precursor[1]) <= mz_distance:
                    found = True
            if not found:
                used_precursors.append(precursor)
        return sum([p[0] for p in used_precursors])

    for protein, precursors in precursors_per_protein_hpp.items():
        score, count = mapping.non_nested_score([(*k,greedy_sequence_precursor_score(v)) for k,v in precursors.items()])
        if score != 0:
            hpp_protein_w_scores.append(fdr.ScoredElement(protein,'XXX_' in protein, score))
            hpp_score_dict[protein] = score

    hpp_fdr_dict = {}
    if len(hpp_protein_w_scores) > 0:
        hpp_fdr_dict = fdr.calculate_fdr(hpp_protein_w_scores, fdr.default_decoy_to_target_function if args.main_fdr == 'picked_hpp' else None)
    
    transform_protein = lambda p: p.replace('XXX_','') if args.leftover_fdr == 'picked_leftover' else p

    for protein, fdr_val in hpp_fdr_dict.items():
        if fdr_val <= 0.01:
            seen_picked.add(transform_protein(protein))

    for protein, precursors in precursors_per_protein_all.items():
        score = sum([greedy_sequence_precursor_score(v) for k,v in precursors.items()])
        #remove corresponding found targets/decoys before running the picked FDR
        if hpp_fdr_dict.get(protein,1) > 0.01 and transform_protein(protein) not in seen_picked:
            if score != 0:
                leftover_protein_w_scores.append(fdr.ScoredElement(protein,'XXX_' in protein, score))
                leftover_score_dict[protein] = score
        if score != 0:
            picked_protein_w_scores.append(fdr.ScoredElement(protein,'XXX_' in protein, score))
            picked_score_dict[protein] = score
            common_protein_w_scores.append(fdr.ScoredElement(protein,'XXX_' in protein, score))
            common_score_dict[protein] = score

    leftover_fdr_dict = {}
    if len(leftover_protein_w_scores) > 0:
        #do proteomicsDB style FDR for canonical FDR
        leftover_fdr_dict = fdr.calculate_fdr(leftover_protein_w_scores, fdr.default_decoy_to_target_function if args.leftover_fdr == 'picked_leftover' else None)

    picked_fdr_dict = {}
    if len(picked_protein_w_scores) > 0:
        #do proteomicsDB style FDR for canonical FDR
        picked_fdr_dict = fdr.calculate_fdr(picked_protein_w_scores, fdr.default_decoy_to_target_function)

    common_fdr_dict = {}
    if len(common_protein_w_scores) > 0:
        #do common FDR as well
        common_fdr_dict = fdr.calculate_fdr(common_protein_w_scores, None)

    if args.output_peptides:
        with open(args.output_peptides,'w') as w:
            header = ['picked_protein_fdr','hpp_protein_fdr','precursor_fdr','psm_fdr','protein_full','protein','protein_type','gene','decoy','all_proteins','pe','ms_evidence','aa_total','database_filename','database_scan','database_usi','sequence','sequence_unmodified','sequence_unmodified_il','charge','score','modifications','pass','type','parent_mass','cosine_filename','cosine_scan','cosine_usi','synthetic_filename','synthetic_scan','synthetic_usi','synthetic_sequence','cosine','synthetic_match','cosine_score_match','explained_intensity','matched_ions','hpp_match','gene_unique','canonical_matches','all_proteins_w_coords','aa_start','aa_end','frag_tol', 'total_unique_exons_covered', 'exons_covered_no_junction', 'exon_junctions_covered', 'all_mapped_exons','datasets','tasks']
            header += ['precursor_fdr_cutoff', 'picked_protein_fdr_cutoff', 'hpp_protein_fdr_cutoff', 'cosine_cutoff', 'explained_intensity_cutoff', 'annotated_ions_cutoff']
            r = csv.DictWriter(w, delimiter = '\t', fieldnames = header, restval='N/A', extrasaction='ignore')
            r.writeheader()
            for precursor in all_precursors:
                proteins = [p for p in precursor['protein'].split(' ###')[0].split(';') if p != '']
                if float(precursor['precursor_fdr']) < 1 and len(proteins) == 1 and ('Canonical' in precursor.get('protein_type','') or 'Contaminant' in precursor.get('protein_type','')):
                    precursor['picked_protein_fdr'] = min(picked_fdr_dict.get(proteins[0],1),1)
                    precursor['hpp_protein_fdr'] = min(hpp_fdr_dict.get(proteins[0],1),1)
                else:
                    precursor['picked_protein_fdr'] = 1
                    precursor['hpp_protein_fdr'] = 1
                precursor.update({
                    'cosine_cutoff':args.cosine_cutoff,
                    'explained_intensity_cutoff':args.explained_intensity_cutoff,
                    'annotated_ions_cutoff':args.annotated_ions_cutoff,
                    'precursor_fdr_cutoff':args.precursor_fdr,
                    'picked_protein_fdr_cutoff':args.picked_protein_fdr,
                    'hpp_protein_fdr_cutoff':args.hpp_protein_fdr
                })
                r.writerow(precursor)

    with open(args.output_mappings, 'w') as w:
        header = ['sequence','protein','protein_full','protein_type','gene','start_aa','end_aa','mismatch_position','protein_aa','peptide_aa','delta_mass','precursor_count','psm_count','best_precursor_usi','best_precursor_filename','best_precursor_scan','best_precursor_charge','best_precursor_sequence']
        r = csv.DictWriter(w, delimiter = '\t', fieldnames = header, restval='N/A', extrasaction='ignore')
        r.writeheader()
        # ProteinMapping = namedtuple('ProteinMapping','protein_accession, start_pos, end_pos, il_ambiguous, mismatches')
        for sequence, mapping_info in pep_mapping_info.items():
            protein_mappings = mapping.string_to_protein_mappings(mapping_info['all_proteins_w_coords'])
            for pmapping in protein_mappings:
                precursor_count = sum([len(position.keys()) for position in frequency[pmapping.protein_accession].values()])
                psm_count = sum([sum(list(position.values())) for position in frequency[pmapping.protein_accession].values()])

                mapped_precursors_at_pos = []
                mapped_psms_at_pos = 0
                best_precursor = None

                if (pmapping.start_pos, pmapping.end_pos) in frequency[pmapping.protein_accession]:
                    for precursor, count in frequency[pmapping.protein_accession][(pmapping.start_pos, pmapping.end_pos)].items():
                        mapped_precursors_at_pos.append(precursor)
                        mapped_psms_at_pos += count

                if len(mapped_precursors_at_pos) > 0 and mapped_precursors_at_pos[0] in representative_per_precursor:
                    best_precursor = representative_per_precursor[mapped_precursors_at_pos[0]]
                    for precursor_sequence_charge in mapped_precursors_at_pos[1:]:
                        if precursor_sequence_charge in representative_per_precursor:
                            precursor = representative_per_precursor[precursor_sequence_charge]
                            if precursor['score'] > best_precursor['score']:
                                best_precursor = precursor

                    delta_mass = 0
                    mismatch_string = ''

                    if len(pmapping.mismatches) > 0:
                        mismatch = pmapping.mismatches[0]
                    else:
                        mismatch = mapping.Mismatch(None,'N/A','N/A')

                    r.writerow({
                        'sequence':sequence,
                        'protein':pmapping.protein_accession,
                        'protein_type':protein_type(pmapping.protein_accession, proteome),
                        'gene':proteome.proteins[pmapping.protein_accession].gene,
                        'start_aa':pmapping.start_pos,
                        'end_aa': pmapping.end_pos,
                        'mismatch_position':int(mismatch.position)+int(pmapping.start_pos) if mismatch.position else 'N/A',
                        'protein_aa':mismatch.protein_aa,
                        'peptide_aa':mismatch.peptide_aa,
                        'delta_mass':unit_delta(mismatch.peptide_aa,mismatch.protein_aa),
                        'precursor_count':len(mapped_precursors_at_pos),
                        'psm_count':mapped_psms_at_pos,
                        'best_precursor_usi':best_precursor['database_usi'] if best_precursor else 'N/A',
                        'best_precursor_filename':best_precursor['database_filename'] if best_precursor else 'N/A',
                        'best_precursor_scan':best_precursor['database_scan'] if best_precursor else 'N/A',
                        'best_precursor_charge':best_precursor['charge'] if best_precursor else 'N/A',
                        'best_precursor_sequence':best_precursor['sequence'] if best_precursor else 'N/A'
                    })

    with open(args.output_proteins, 'w') as fo:

        overlap_keys = find_overlap({},{},0,0,"")[0].keys()

        fieldnames = [
            'protein',
            'protein_full',
            'protein_type',
            'gene',
            'pe',
            'aa_total',
            'hpp_score',
            'hpp_fdr',
            'picked_score',
            'picked_fdr',
            'common_score',
            'common_fdr',
            'leftover_score',
            'leftover_fdr',
            'num_sequences_incl_shared',
            'num_sequences',
            'coverage_incl_shared',
            'precursor_fdr_cutoff',
            'hpp_protein_fdr_cutoff',
            'picked_protein_fdr_cutoff',
            'cosine_cutoff',
            'explained_intensity_cutoff',
            'annotated_ions_cutoff'
        ] + ['_dyn_#neXtProt Release {}'.format(release) for release in sorted(nextprot_releases_pe.keys())]

        for overlap_suffix in ['','_w_synthetic','_w_synthetic_cosine','_just_current','_just_current_iso_unique']:
            fieldnames.extend(['{}{}'.format(k,overlap_suffix) for k in overlap_keys])

        w = csv.DictWriter(fo, delimiter = '\t', fieldnames = fieldnames, extrasaction='ignore')
        w.writeheader()

        hpp_per_protein = {}
        coverage_per_protein = {}

        for i,(protein,protein_entry) in enumerate(proteome.proteins.items()):

            is_canonical = protein_entry.db == 'sp' and not protein_entry.iso

            protein_dict = {
                'protein': protein,
                'protein_full':mapping.protein_to_full_uniprot_string(proteome.proteins[protein]),
                'protein_type':protein_type(protein, proteome),
                'pe': nextprot_pe.get(protein_entry.id,0) if is_canonical else 0,
                'aa_total':protein_entry.length,
                'gene':protein_entry.gene,
                'hpp_score':hpp_score_dict.get(protein,0),
                'hpp_fdr':min(1,hpp_fdr_dict.get(protein,1)),
                'picked_score':picked_score_dict.get(protein,0),
                'picked_fdr':min(1,picked_fdr_dict.get(protein,1)),
                'common_score':common_score_dict.get(protein,0),
                'common_fdr':min(1,common_fdr_dict.get(protein,1)),
                'leftover_score':leftover_score_dict.get(protein,0),
                'leftover_fdr':min(1,leftover_fdr_dict.get(protein,1)),
                'num_sequences':len(precursors_per_protein_all.get(protein,[])),
                'num_sequences_incl_shared':len(precursors_per_protein_non_unique.get(protein,[])),
            }

            pass_comparison_picked_fdr,pass_comparison_hpp_fdr = comparison_picked_fdr.get(protein,1) <= args.picked_protein_fdr_comparison, comparison_hpp_fdr.get(protein,1) <= args.hpp_protein_fdr_comparison 
            pass_picked_fdr, pass_hpp_fdr = protein_dict['picked_fdr'] <= args.picked_protein_fdr, protein_dict['hpp_fdr'] <= args.hpp_protein_fdr 
            if args.main_fdr == 'traditional':
                pass_picked_fdr, pass_hpp_fdr = protein_dict['common_fdr'] <= args.picked_protein_fdr, protein_dict['common_fdr'] <= args.hpp_protein_fdr 

            if not args.filter_rows_fdr or any([pass_comparison_picked_fdr, pass_comparison_hpp_fdr, pass_picked_fdr, pass_hpp_fdr]):

                all_protein_mappings = protein_mapping[protein]

                compare_mappings = defaultdict(dict)

                for sequence, mappings in all_protein_mappings.items():
                    found = sequences_found.get(sequence)
                    if found:
                        if found.added.match:
                            compare_mappings['all'].update({sequence:mappings})
                        if found.added.synthetic_match:
                            compare_mappings['all_synthetic_match'].update({sequence:mappings})
                        if found.added.synthetic_match_cosine:
                            compare_mappings['all_synthetic_match_cosine'].update({sequence:mappings})
                        if found.hpp:
                            if found.comparison.match:
                                compare_mappings['comparison_match'].update({sequence:mappings})
                            if found.comparison.synthetic_match:
                                compare_mappings['comparison_synthetic_match'].update({sequence:mappings})
                            if found.comparison.synthetic_match_cosine:
                                compare_mappings['comparison_synthetic_match_cosine'].update({sequence:mappings})
                            if found.added.match:
                                compare_mappings['added_match'].update({sequence:mappings})
                            if found.added.synthetic_match:
                                compare_mappings['added_synthetic_match'].update({sequence:mappings})
                            if found.added.synthetic_match_cosine:
                                compare_mappings['added_synthetic_match_cosine'].update({sequence:mappings})
                            if found.isoform_unique:
                                compare_mappings['isoform_unique'].update({sequence:mappings})

                protein_dict['coverage_incl_shared'] = mapping.find_coverage([],list(compare_mappings['all'].values()), protein_entry.length)/protein_entry.length

                if is_canonical:
                    protein_dict.update(find_overlap(compare_mappings['comparison_match'],compare_mappings['added_match'],int(protein_dict['aa_total']),int(protein_dict['pe']),'',10,pass_picked_fdr, pass_hpp_fdr, pass_comparison_picked_fdr,pass_comparison_hpp_fdr)[0])
                else:
                    protein_dict.update(find_overlap(compare_mappings['comparison_match'],compare_mappings['isoform_unique'],int(protein_dict['aa_total']),int(protein_dict['pe']),'',10,pass_picked_fdr, pass_hpp_fdr, pass_comparison_picked_fdr,pass_comparison_hpp_fdr)[0])

                #comparison_proteins.get(protein,{})
                protein_dict.update(find_overlap(compare_mappings['comparison_synthetic_match'],compare_mappings['added_synthetic_match'],int(protein_dict['aa_total']),int(protein_dict['pe']),'_w_synthetic',10,pass_picked_fdr, pass_hpp_fdr, pass_comparison_picked_fdr,pass_comparison_hpp_fdr)[0])
                protein_dict.update(find_overlap(compare_mappings['comparison_synthetic_match_cosine'],compare_mappings['added_synthetic_match_cosine'],int(protein_dict['aa_total']),int(protein_dict['pe']),'_w_synthetic_cosine',10,pass_picked_fdr, pass_hpp_fdr, pass_comparison_picked_fdr,pass_comparison_hpp_fdr)[0])
                if is_canonical:
                    protein_dict.update(find_overlap({},compare_mappings['added_match'],int(protein_dict['aa_total']),int(protein_dict['pe']),'_just_current',10,pass_picked_fdr, pass_hpp_fdr)[0])
                else:
                    protein_dict.update(find_overlap({},{},int(protein_dict['aa_total']),int(protein_dict['pe']),'_just_current',10,pass_picked_fdr, pass_hpp_fdr)[0])

                hpp_per_protein[protein] = protein_dict['combined_hpp_just_current']
                coverage_per_protein[protein] = protein_dict['coverage_incl_shared']

                protein_dict.update(find_overlap({},compare_mappings['isoform_unique'],int(protein_dict['aa_total']),int(protein_dict['pe']),'_just_current_iso_unique',10,pass_picked_fdr, pass_hpp_fdr)[0])

                protein_dict.update({
                    'cosine_cutoff':args.cosine_cutoff,
                    'explained_intensity_cutoff':args.explained_intensity_cutoff,
                    'annotated_ions_cutoff':args.annotated_ions_cutoff,
                    'precursor_fdr_cutoff':args.precursor_fdr,
                    'picked_protein_fdr_cutoff':args.picked_protein_fdr,
                    'hpp_protein_fdr_cutoff':args.hpp_protein_fdr
                })

                for release, pe_dict in nextprot_releases_pe.items():
                    protein_dict['_dyn_#neXtProt Release {}'.format(release)] = pe_dict.get(protein_entry.id,0) if is_canonical else 0
                w.writerow(protein_dict)

    with open(args.output_dataset_proteins_hpp, 'w') as f_dataset_hpp, open(args.output_dataset_proteins_all, 'w') as f_dataset_all, open(args.output_task_proteins_hpp, 'w') as f_task_hpp, open(args.output_task_proteins_all, 'w') as f_task_all:

        dataset_header = lambda d: '_dyn_#{}'.format("{}/{}".format(pxd_mapping[d],d) if d in pxd_mapping else d)

        def get_task_name(task_id):
            url = 'http://proteomics2.ucsd.edu/ProteoSAFe/status_json.jsp?task={}'.format(task_id)
            return requests.get(url, verify=False).json()['description'].split(' [launched')[0]

        task_names = {}

        for task in all_tasks:
            task_names[task] = get_task_name(task)

        task_header = lambda t: '_dyn_#{}'.format(task_names[t])

        common_fieldnames = [
            'protein',
            'protein_full',
            'protein_type',
            'gene',
            'pe',
            'aa_total',
            'hpp_sequences',
            'hpp_score',
            'hpp_fdr',
            'picked_score',
            'picked_fdr',
            'common_score',
            'common_fdr',
            'leftover_score',
            'leftover_fdr',
            'coverage_incl_shared',
        ] + ['_dyn_#neXtProt Release {}'.format(release) for release in sorted(nextprot_releases_pe.keys())]
        
        dataset_fieldnames = common_fieldnames + [dataset_header(d) for d in all_datasets]
        task_fieldnames = common_fieldnames + [task_header(t) for t in all_tasks]

        w_dataset_hpp = csv.DictWriter(f_dataset_hpp, delimiter = '\t', fieldnames = dataset_fieldnames, restval=0, extrasaction='ignore')
        w_dataset_hpp.writeheader()

        w_dataset_all = csv.DictWriter(f_dataset_all, delimiter = '\t', fieldnames = dataset_fieldnames, restval=0, extrasaction='ignore')
        w_dataset_all.writeheader()

        w_task_hpp = csv.DictWriter(f_task_hpp, delimiter = '\t', fieldnames = task_fieldnames, restval=0, extrasaction='ignore')
        w_task_hpp.writeheader()

        w_task_all = csv.DictWriter(f_task_all, delimiter = '\t', fieldnames = task_fieldnames, restval=0, extrasaction='ignore')
        w_task_all.writeheader()

        for i,(protein,protein_entry) in enumerate(proteome.proteins.items()):

            if protein in hpp_per_protein:

                is_canonical = protein_entry.db == 'sp' and not protein_entry.iso

                protein_dict = {
                    'protein': protein,
                    'protein_full':mapping.protein_to_full_uniprot_string(proteome.proteins[protein]),
                    'protein_type':protein_type(protein, proteome),
                    'pe': nextprot_pe.get(protein_entry.id,0) if is_canonical else 0,
                    'aa_total':protein_entry.length,
                    'hpp_sequences': hpp_per_protein.get(protein,0),
                    'gene':protein_entry.gene,
                    'hpp_score':hpp_score_dict.get(protein,0),
                    'hpp_fdr':min(1,hpp_fdr_dict.get(protein,1)),
                    'picked_score':picked_score_dict.get(protein,0),
                    'picked_fdr':min(1,picked_fdr_dict.get(protein,1)),
                    'common_score':common_score_dict.get(protein,0),
                    'common_fdr':min(1,common_fdr_dict.get(protein,1)),
                    'leftover_score':leftover_score_dict.get(protein,0),
                    'leftover_fdr':min(1,leftover_fdr_dict.get(protein,1)),
                    'coverage_incl_shared':coverage_per_protein.get(protein,0)
                }

                pass_picked_fdr, pass_hpp_fdr = protein_dict['picked_fdr'] <= args.picked_protein_fdr, protein_dict['hpp_fdr'] <= args.hpp_protein_fdr 

                all_protein_mappings = protein_mapping[protein]

                sequences_per_dataset_all = defaultdict(dict)
                sequences_per_task_all = defaultdict(dict)
                sequences_per_dataset_hpp = defaultdict(dict)
                sequences_per_task_hpp = defaultdict(dict)

                for sequence, mappings in all_protein_mappings.items():
                    found = sequences_found.get(sequence)
                    if found:
                        if found.added.match:
                            for dataset in datasets_per_sequence.get(sequence,set()):
                                sequences_per_dataset_all[dataset].update({sequence:mappings})
                                if found.hpp:
                                    sequences_per_dataset_hpp[dataset].update({sequence:mappings})
                            for task in tasks_per_sequence.get(sequence,set()):
                                sequences_per_task_all[task].update({sequence:mappings})
                                if found.hpp:
                                    sequences_per_task_hpp[task].update({sequence:mappings})

                for release, pe_dict in nextprot_releases_pe.items():
                    protein_dict['_dyn_#neXtProt Release {}'.format(release)] = pe_dict.get(protein_entry.id,0) if is_canonical else 0

                protein_dict_dataset_hpp = protein_dict.copy()
                protein_dict_dataset_all = protein_dict.copy()

                protein_dict_task_hpp = protein_dict.copy()
                protein_dict_task_all = protein_dict.copy()

                for dataset,positions in sequences_per_dataset_hpp.items():
                    overlaps = find_overlap({},positions,int(protein_dict['aa_total']),int(protein_dict['pe']),'', 10, True, True, False, False)[0]
                    protein_dict_dataset_hpp[dataset_header(dataset)] = overlaps['combined_hpp']

                for dataset,positions in sequences_per_dataset_all.items():
                    protein_dict_dataset_all[dataset_header(dataset)] = len(positions)

                for task,positions in sequences_per_task_hpp.items():
                    overlaps = find_overlap({},positions,int(protein_dict['aa_total']),int(protein_dict['pe']),'', 10, True, True, False, False)[0]
                    protein_dict_task_hpp[task_header(task)] = overlaps['combined_hpp']

                for task,positions in sequences_per_task_all.items():
                    protein_dict_task_all[task_header(task)] = len(positions)

                w_dataset_hpp.writerow(protein_dict_dataset_hpp)
                w_dataset_all.writerow(protein_dict_dataset_all)
                w_task_hpp.writerow(protein_dict_task_hpp)
                w_task_all.writerow(protein_dict_task_all)

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

        w = csv.DictWriter(fo, delimiter = '\t', fieldnames = fieldnames, extrasaction='ignore')
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
        for protein, peptides in protein_mapping.items():
            if 'XXX_' not in protein:
                protein_dict = {
                    'protein': protein,
                    'pe': nextprot_pe.get(protein,0),
                    'aa_total':proteome.proteins.get(protein).length,
                    'hpp_score':hpp_score_dict.get(protein,0),
                    'hpp_fdr':min(1,hpp_fdr_dict.get(protein,1)),
                    'picked_score':picked_score_dict.get(protein,0),
                    'picked_fdr':min(1,picked_fdr_dict.get(protein,1))
                }
                pass_picked_fdr, pass_hpp_fdr = protein_dict['picked_fdr'] <= args.picked_protein_fdr, protein_dict['hpp_fdr'] <= args.hpp_protein_fdr 
                hupo_mapping[protein] = find_overlap({},peptides,int(protein_dict['aa_total']),int(protein_dict['pe']),'', 10, pass_picked_fdr, pass_hpp_fdr,0,0)[0]['new_hpp']
        for protein, peptides in protein_mapping.items():
            if 'XXX_' not in protein:
                unique_mapping[protein] = len(set([peptide for peptide, mappings in peptides.items()]))
        explorer_export.output_for_explorer(args.explorers_output, pep_mapping_info, representative_per_precursor, args.external_provenance, hupo_mapping, unique_mapping, args.library_name, args.library_version)

if __name__ == '__main__':
    main()
