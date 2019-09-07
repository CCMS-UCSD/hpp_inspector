#!/usr/bin/python

"""

PSM Utilities to read psms

"""

import ming_fileio_library
import math
import re

class PSM:
    def __init__(self, filename, scan, annotation, score, decoy, protein):
        self.filename = filename
        self.scan = scan
        self.annotation = annotation
        self.score = score
        self.decoy = decoy
        self.protein = protein
        self.fdr = -1.0

    def __str__(self):
        return self.filename + ":" + self.scan + " " + self.annotation

    def __repr__(self):
        return str(self)

    def is_decoy(self):
        return self.decoy

    def sorting_value(self):
        return self.score

    def get_stripped_sequence(self):
        sequence = self.annotation
        p = re.compile('\W|\d')
        sequence = p.sub("", sequence)
        return sequence


class PSMset:
    def __init__(self, name):
        self.name = name
        self.psms = []

    def __len__(self):
        return len(self.psms)

    #Loading a TSV File from MSGFDB
    def load_MSGF_tsvfile(self, filename):
        self.psms += parse_MSGF_tsvfile(filename)

    #Filter PSMs to given FDR
    def filter_to_fdr(self, fdr):
        filtered_psms = filter_psm_fdr(self.psms, fdr)
        print "Filtered " + str(len(self.psms)) + " to " + str(len(filtered_psms))
        self.psms = filtered_psms

    def filter_to_fdr_by_length(self, fdr):
        output_psms = []
        peptide_length_map = {}
        for psm in self.psms:
            peptide_length = len(psm.get_stripped_sequence())
            if not peptide_length in peptide_length_map:
                peptide_length_map[peptide_length] = []
            peptide_length_map[peptide_length].append(psm)

        for peptide_length in peptide_length_map:
            filtered_psms = filter_psm_fdr(peptide_length_map[peptide_length], fdr)
            print "Filtered Length " + str(peptide_length) + " " + str(len(peptide_length_map[peptide_length])) + " to " + str(len(filtered_psms))
            output_psms += filtered_psms
        self.psms = output_psms

    #Calculate FDR of PSM Set
    def calculate_fdr(self):
        running_target_count = 0
        running_decoy_count = 0
        for psm in self.psms:
            if psm.is_decoy() == 0:
                running_target_count += 1
            else:
                running_decoy_count += 1

        current_fdr = float(running_decoy_count) / float(running_target_count)
        return current_fdr


class PeptideVariant:
    def __init__(self, variant_sequence):
        self.variant_sequence = variant_sequence
        self.psms = []

    def __str__(self):
        return self.variant_sequence + "\t" + str(self.sorting_value()) + "\t" + str(self.is_decoy())

    def add_psm(self, psm_object):
        self.psms.append(psm_object)

    def is_decoy(self):
        return self.psms[0].is_decoy()

    def sorting_value(self):
        max_score = 0
        for psm in self.psms:
            max_score = max(psm.sorting_value(), max_score)
        return max_score

    def get_spectrum_count(self):
        return len(self.psms)

    def get_stripped_sequence(self):
        sequence = self.variant_sequence
        p = re.compile('\W|\d')
        sequence = p.sub("", sequence)
        return sequence

    def sequence_length(self):
        return len(self.get_stripped_sequence())

    #Research-y Portion

###
 # Class to hold a set of library peptides
###
class PeptideVariantSet:
    def __init__(self, name):
        self.name = name
        self.peptide_list = []
        self.peptide_map = {}

    def __len__(self):
        return len(self.peptide_list)

    def get_total_spectra_count(self):
        return sum(self.get_spectra_count_list())
        total_ms_ms_count = 0
        for variant in self.peptide_list:
            total_ms_ms_count += len(variant.psms)
        return total_ms_ms_count

    #Get total peptides regardless of modifications
    def get_total_unique_sequence_count(self):
        return len(self.get_unique_sequences())

    def get_unique_sequences_spectrum_count_map(self):
        sequence_map = {}
        for variant in self.peptide_list:
            sequence = variant.variant_sequence
            p = re.compile('\W|\d')
            sequence = p.sub("", sequence)
            if not sequence in sequence_map:
                sequence_map[sequence] = 0
            sequence_map[sequence] += variant.get_spectrum_count()

        return sequence_map

    def get_unique_sequences(self):
        sequence_map = {}
        for variant in self.peptide_list:
            sequence = variant.variant_sequence
            p = re.compile('\W|\d')
            sequence = p.sub("", sequence)
            sequence_map[sequence] = 1

        return sequence_map.keys()

    #Returns a list of spectral counts for each variant
    def get_spectra_count_list(self):
        ms_ms_count_list = []
        for variant in self.peptide_list:
            ms_ms_count_list.append(len(variant.psms))
        return ms_ms_count_list

    def add_psms_set(self, psm_set):
        self.add_psms_list(psm_set.psms)

    def add_psms_list(self, psm_list):
        for psm in psm_list:
            if not(psm.annotation in self.peptide_map):
                peptide_variant = PeptideVariant(psm.annotation)
                self.peptide_list.append(peptide_variant)
                self.peptide_map[psm.annotation] = peptide_variant
            self.peptide_map[psm.annotation].add_psm(psm)

    #Appending a variant set
    def add_variant_set(self, variant_set):
        for variant in variant_set.peptide_list:
            if not(variant.variant_sequence in self.peptide_map):
                self.peptide_list.append(variant)
                self.peptide_map[variant.variant_sequence] = variant
            else:
                self.add_psms_list(variant.psms)

    def add_variant(self, variant_obj):
        if not(variant_obj.variant_sequence in self.peptide_map):
            self.peptide_list.append(variant_obj)
            self.peptide_map[variant_obj.variant_sequence] = variant_obj
        else:
            self.add_psms_list(variant_obj.psms)

    def remove_variant(self, variant_obj):
        self.peptide_list.remove(variant_obj)
        del self.peptide_map[variant_obj.variant_sequence]

    def filter_to_fdr(self, fdr):
        filtered_peptides = filter_psm_fdr(self.peptide_list, fdr)
        print "Filtered " + str(len(self.peptide_list)) + " to " + str(len(filtered_peptides))
        self.peptide_list = filtered_peptides
        self.peptide_map = {}
        for variant in self.peptide_list:
            self.peptide_map[variant.variant_sequence] = variant

    def calculate_fdr(self):
        running_target_count = 0
        running_decoy_count = 0
        for psm in self.peptide_list:
            if psm.is_decoy() == 0:
                running_target_count += 1
            else:
                running_decoy_count += 1

        current_fdr = float(running_decoy_count) / float(running_target_count)
        return current_fdr




###
 # Takes as input a filename
 # Returns a list of PSM
###
def parse_MSGF_tsvfile(filename):
    rows, table_data = ming_fileio_library.parse_table_with_headers(filename)
    scan_header = "Scan#"
    peptide_header = "Peptide"
    protein_header = "Protein"
    score_header = "P-value"
    filename_header = "#SpecFile"
    charge_header = "Charge"

    decoy_indicator = "REV_"

    psm_list = []

    for i in range(rows):
        scan = table_data[scan_header][i]
        peptide = table_data[peptide_header][i]
        protein = table_data[protein_header][i]
        score = -math.log10(float(table_data[score_header][i]))
        #print table_data[score_header][i] + "\t" + str(score)
        filename = table_data[filename_header][i]
        charge = table_data[charge_header][i]
        decoy = 0

        #Stripping peptide dots
        if peptide[1] == "." and peptide[-2] == ".":
            peptide = peptide[2:-2]


        if protein.find(decoy_indicator) != -1:
            decoy = 1

        #Adding charge state to peptide name
        peptide += "." + str(charge)

        new_psm = PSM(filename, scan, peptide, score, decoy, protein)
        psm_list.append(new_psm)

    return psm_list


###
 # Filtering PSM results so that the returned set is at the
 # prescribed FDR
###
def filter_psm_fdr(input_psms, fdr_percentage):
    #TODO Add Sorting
    #Sorting
    #print "Sorting shit to " + str(fdr_percentage*100) + "%"
    #input_psms = sorted(input_psms, key=PSM.sorting_value)
    input_psms = sorted(input_psms, key=lambda psm: psm.sorting_value(), reverse=True)

    running_target_count = 1
    running_decoy_count = 0

    output_psms = []

    for psm in input_psms:
        if psm.is_decoy() == 0:
            running_target_count += 1
        else:
            running_decoy_count += 1

        current_fdr = float(running_decoy_count) / float(running_target_count)

        psm.fdr = current_fdr

    #Properly finding the min q value for PSM
    min_fdr = 1
    for i in range(len(input_psms)):
        index_to_check = len(input_psms) - i - 1
        min_fdr = min(min_fdr, input_psms[index_to_check].fdr)
        input_psms[index_to_check].fdr = min_fdr

    for psm in input_psms:
        if psm.fdr < fdr_percentage:
            output_psms.append(psm)

    return output_psms
