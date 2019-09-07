#!/usr/bin/python


import sys
import getopt
import os
import ming_protein_library
import copy


if (sys.version_info > (3, 0)):
    # Python 3 code in this block
    import _pickle as pickle
else:
    # Python 2 code in this block
    import cPickle as pickle




def load_peptide_list(input_peptide_filename):
    all_peptides = []
    for line in open(input_peptide_filename):
        peptide = line.rstrip()
        all_peptides.append(peptide)
    return all_peptides

#Looking for modifications, only looking at point annotations
def modification_residue_statistics(modification_list, coverage_array):
    modifications_covered = 0
    modifications_total = 0
    for annotation in modification_list:
        if "position" in annotation["location"][0].keys():
            try:
                residue_number = int(annotation["location"][0]["position"][0]["@position"])
                residue_index = residue_number - 1
                if residue_index >= len(coverage_array):
                    continue
                if coverage_array[residue_index] == 1:
                    modifications_covered += 1
                modifications_total += 1
            except KeyboardInterrupt:
                raise
            except:
                continue

    if modifications_total == 0:
        return 0,0,0

    percentage_covered = float(modifications_covered)/float(modifications_total)
    return modifications_covered, modifications_total, percentage_covered

#Looking for disulfide bond sites, specifies locations of the two cysteines that form the cystine
def disulfide_bond_residue_statistics(modification_list, coverage_array):
    disulfide_bond_sites = []
    for annotation in modification_list:
        try:
            if "position" in annotation["location"][0].keys():
                position = int(annotation["location"][0]["position"][0]["@position"])
                disulfide_bond_sites.append(position - 1)

            if "begin" in annotation["location"][0].keys():
                end_site = int(annotation["location"][0]["end"][0]["@position"])
                start_site = int(annotation["location"][0]["begin"][0]["@position"])
                disulfide_bond_sites.append(end_site - 1)
                disulfide_bond_sites.append(start_site - 1)
        except KeyboardInterrupt:
            raise
        except:
            continue



    disulfide_bond_sites = list(set(disulfide_bond_sites))

    covered_site_count = 0
    for site_index in disulfide_bond_sites:
        if site_index >= len(coverage_array):
            continue
        if coverage_array[site_index] == 1:
            covered_site_count += 1

    percentage_covered = float(covered_site_count)/float(len(disulfide_bond_sites))

    return covered_site_count, len(disulfide_bond_sites), percentage_covered

#Looking at the SNPs and seeing the proteins that result from that with one natural variant at a time
def natural_variant_statistics(modification_list, protein_object, peptide_list):
    site_total = 0
    covered_site_count = 0
    #Doing stuff
    for annotation in modification_list:
        if "position" in annotation["location"][0].keys():
            try:
                position = int(annotation["location"][0]["position"][0]["@position"])
                residue_number = position
                residue_index = residue_number - 1

                if not "original" in annotation.keys():
                    continue

                substitution_from = annotation["original"][0]
                substitution_to = annotation["variation"][0]

                #print substitution_from + "->" + substitution_to

                #Creates variant of protein sequence
                sequence_list = list(protein_object.sequence)
                sequence_list[residue_index] = substitution_to
                updated_protein_sequence = "".join(sequence_list)

                #print updated_protein_sequence

                site_total += 1

                #Create new protein object
                snp_protein = ming_protein_library.Protein("DB Coverage", updated_protein_sequence, "GENE")
                snp_coverage = snp_protein.coverage_by_amino_acids(peptide_list)
                if residue_index >= len(coverage_array):
                    continue
                if snp_coverage[residue_index] == 1:
                    covered_site_count += 1
            except KeyboardInterrupt:
                raise
            except:
                continue



    if site_total == 0:
        return 0,0,0

    percentage_covered = float(covered_site_count)/float(site_total)

    return covered_site_count, site_total, percentage_covered

#Isoforms as listed by uniprot, todo, handle multiple changes to isoforms, not totally correct at the moment
def isoforms_statistics(modification_list, protein_object):
    for annotation in modification_list:
        #Getting the note
        note_string = annotation.annotation_text.split(";")[1]

        note_value = note_string.split("=")[1]

        #Getting the isoform number
        isoform_number = note_value.split(".")[0][-1]

        #Isoform from
        if note_value.split(".")[1].strip() == "Missing":
            continue

        from_sequence = note_value.split(".")[1].strip().split("->")[0]
        to_sequence = note_value.split(".")[1].strip().split("->")[1]


        substitution_start_index = annotation.start_site - 1
        substitution_end_index = annotation.end_site - 1

        if from_sequence != protein_object.sequence[substitution_start_index: substitution_end_index + 1]:
            continue


        print( "ISOFORM: " + str(isoform_number))
        print(from_sequence + "->" + to_sequence)


    return 0, 0, 0


#Looking for number of junction peptides
def calculate_protein_statistics(uniprot_accession, protein_to_annotations_map, protein_object, peptide_list, coverage_array):
    return_map = {}
    return_map["mod_covered"] = 0
    return_map["mod_total"] = 0
    return_map["mod_percentage"] = 0

    return_map["sites_covered"] = 0
    return_map["sites_total"] = 0
    return_map["sites_percentage"] = 0

    return_map["activesites_covered"] = 0
    return_map["activesites_total"] = 0
    return_map["activesites_percentage"] = 0

    return_map["bindingsites_covered"] = 0
    return_map["bindingsites_total"] = 0
    return_map["bindingsites_percentage"] = 0

    return_map["disulfide_covered"] = 0
    return_map["disulfide_total"] = 0
    return_map["disulfide_percentage"] = 0

    return_map["glycosylation_covered"] = 0
    return_map["glycosylation_total"] = 0
    return_map["glycosylation_percentage"] = 0

    return_map["metalbinding_covered"] = 0
    return_map["metalbinding_total"] = 0
    return_map["metalbinding_percentage"] = 0

    return_map["snp_covered"] = 0
    return_map["snp_total"] = 0
    return_map["snp_percentage"] = 0

    return_map["protein_name"] = "Protein Name"
    return_map["gene"] = "gene"

    if not uniprot_accession in protein_to_annotations_map:
        return return_map


    try:
        if "name" in protein_to_annotations_map[uniprot_accession]:
            return_map["protein_name"] = str(protein_to_annotations_map[uniprot_accession]["name"])

        if "gene" in protein_to_annotations_map[uniprot_accession]:
            return_map["gene"] = str(protein_to_annotations_map[uniprot_accession]["gene"])

        protein_annotation_map = protein_to_annotations_map[uniprot_accession]["features"]
        if "modified residue" in protein_annotation_map.keys():
            modifications_covered, modifications_total, percentage_covered = modification_residue_statistics(protein_annotation_map["modified residue"], coverage_array)
            return_map["mod_covered"] = modifications_covered
            return_map["mod_total"] = modifications_total
            return_map["mod_percentage"] = percentage_covered

        if "site" in protein_annotation_map.keys():
            annotations_covered, annotations_total, percentage_covered = modification_residue_statistics(protein_annotation_map["site"], coverage_array)
            return_map["sites_covered"] = annotations_covered
            return_map["sites_total"] = annotations_total
            return_map["sites_percentage"] = percentage_covered

        if "active site" in protein_annotation_map.keys():
            annotations_covered, annotations_total, percentage_covered = modification_residue_statistics(protein_annotation_map["active site"], coverage_array)
            return_map["activesites_covered"] = annotations_covered
            return_map["activesites_total"] = annotations_total
            return_map["activesites_percentage"] = percentage_covered

        if "binding site" in protein_annotation_map.keys():
            annotations_covered, annotations_total, percentage_covered = modification_residue_statistics(protein_annotation_map["binding site"], coverage_array)
            return_map["bindingsites_covered"] = annotations_covered
            return_map["bindingsites_total"] = annotations_total
            return_map["bindingsites_percentage"] = percentage_covered

        if "disulfide bond" in protein_annotation_map.keys():
            annotations_covered, annotations_total, percentage_covered = disulfide_bond_residue_statistics(protein_annotation_map["disulfide bond"], coverage_array)
            return_map["disulfide_covered"] = annotations_covered
            return_map["disulfide_total"] = annotations_total
            return_map["disulfide_percentage"] = percentage_covered

        if "glycosylation site" in protein_annotation_map.keys():
            annotations_covered, annotations_total, percentage_covered = modification_residue_statistics(protein_annotation_map["glycosylation site"], coverage_array)
            return_map["glycosylation_covered"] = annotations_covered
            return_map["glycosylation_total"] = annotations_total
            return_map["glycosylation_percentage"] = percentage_covered

        if "metal ion-binding site" in protein_annotation_map.keys():
            annotations_covered, annotations_total, percentage_covered = modification_residue_statistics(protein_annotation_map["metal ion-binding site"], coverage_array)
            return_map["metalbinding_covered"] = annotations_covered
            return_map["metalbinding_total"] = annotations_total
            return_map["metalbinding_percentage"] = percentage_covered

        if "sequence variant" in protein_annotation_map.keys():
            annotations_covered, annotations_total, percentage_covered = natural_variant_statistics(protein_annotation_map["sequence variant"], protein_object, peptide_list)
            return_map["snp_covered"] = annotations_covered
            return_map["snp_total"] = annotations_total
            return_map["snp_percentage"] = percentage_covered

        #if "Alternative sequence" in protein_annotation_map:
        #    annotations_covered, annotations_total, percentage_covered = isoforms_statistics(protein_annotation_map["Alternative sequence"], protein_object )



    except (KeyboardInterrupt, SystemExit):
        raise
    except:
        #import json
        #print(json.dumps(protein_to_annotations_map[uniprot_accession]))
        raise

        return return_map

    #print return_map

    return return_map

def usage():
    print("<input peptides list> <input uniprot fasta file> <input pkl of gff file> <output stats>")

def main():
    input_peptide_list_filename = sys.argv[1]
    input_uniprot_fasta_filename = sys.argv[2]
    input_uniprot_gff_pkl_filename = sys.argv[3]
    input_nextprot_pe = sys.argv[4]
    output_file = open(sys.argv[5], "w")

    proteome = ming_protein_library.parse_fasta_proteome_file(input_uniprot_fasta_filename,input_nextprot_pe)
    protein_to_annotations_map = pickle.load(open(input_uniprot_gff_pkl_filename, "r"))
    peptide_list = load_peptide_list(input_peptide_list_filename)

    output_header = "Protein\tProteinName\tModCovered\tModTotal\tModPercent\tSiteCovered\tSiteTotal\tSitePercent\t"
    output_header += "ActiveSiteCovered\tActiveSiteTotal\tActiveSitePercent\tBindingSiteCovered\tBindingSiteTotal\tBindingSitePercent\t"
    output_header += "DisulfideSiteCovered\tDisulfideSiteTotal\tDisulfideSitePercent\t"
    output_header += "GlycosylationCovered\tGlycosylationTotal\tGlycosylationPercent\t"
    output_header += "MetalBindingCovered\tMetalBindingTotal\tMetalBindingPercent\t"
    output_header += "SNPCovered\tSNPTotal\tSNPPercent"

    output_file.write(output_header + "\n")

    for protein in proteome.protein_list:
        uniprot_accession = protein.protein.split("|")[1]
        #Skipping alternative isoforms of protein, becasue annotations are incorrect for them
        if uniprot_accession.find("-") != -1:
            continue

        print (uniprot_accession)

        coverage_array = protein.coverage_by_amino_acids(peptide_list)
        protein_stats_map = calculate_protein_statistics(uniprot_accession, protein_to_annotations_map, protein, peptide_list, coverage_array)

        output_string = uniprot_accession + "\t"
        output_string += protein_stats_map["protein_name"] + "\t"

        output_string += str(protein_stats_map["mod_covered"]) + "\t"
        output_string += str(protein_stats_map["mod_total"]) + "\t"
        output_string += str(protein_stats_map["mod_percentage"]) + "\t"

        output_string += str(protein_stats_map["sites_covered"]) + "\t"
        output_string += str(protein_stats_map["sites_total"]) + "\t"
        output_string += str(protein_stats_map["sites_percentage"]) + "\t"

        output_string += str(protein_stats_map["activesites_covered"]) + "\t"
        output_string += str(protein_stats_map["activesites_total"]) + "\t"
        output_string += str(protein_stats_map["activesites_percentage"]) + "\t"

        output_string += str(protein_stats_map["bindingsites_covered"]) + "\t"
        output_string += str(protein_stats_map["bindingsites_total"]) + "\t"
        output_string += str(protein_stats_map["bindingsites_percentage"]) + "\t"

        output_string += str(protein_stats_map["disulfide_covered"]) + "\t"
        output_string += str(protein_stats_map["disulfide_total"]) + "\t"
        output_string += str(protein_stats_map["disulfide_percentage"]) + "\t"

        output_string += str(protein_stats_map["glycosylation_covered"]) + "\t"
        output_string += str(protein_stats_map["glycosylation_total"]) + "\t"
        output_string += str(protein_stats_map["glycosylation_percentage"]) + "\t"

        output_string += str(protein_stats_map["metalbinding_covered"]) + "\t"
        output_string += str(protein_stats_map["metalbinding_total"]) + "\t"
        output_string += str(protein_stats_map["metalbinding_percentage"]) + "\t"

        output_string += str(protein_stats_map["snp_covered"]) + "\t"
        output_string += str(protein_stats_map["snp_total"]) + "\t"
        output_string += str(protein_stats_map["snp_percentage"]) + "\t"

        output_file.write(output_string + "\n")
        #print output_string







if __name__ == "__main__":
    main()
