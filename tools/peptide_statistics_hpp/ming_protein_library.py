#!/usr/bin/python


import sys
import getopt
import os
import re
import ming_fileio_library
import distance
import regex
from collections import defaultdict, namedtuple
from csv import DictReader

if (sys.version_info > (3, 0)):
    # Python 3 code in this block
    import _pickle as pickle
else:
    # Python 2 code in this block
    import cPickle as pickle

ProteomeMatch = namedtuple('ProteomeMatch','protein, start_pos, end_pos, il_ambiguous, mismatches')

class ProteinAnnotation:
    def __init__(self, annotation_class, start_site, end_site, annotation_text):
        self.annotation_class = annotation_class
        self.start_site = start_site
        self.end_site = end_site
        self.annotation_text = annotation_text


def parse_fasta_proteome_file(filename, output_raw = False, external_PE = False):

    pe_ref = {}

    if (external_PE):
        print(external_PE)
        with open(external_PE, 'r') as f:
            r = DictReader(f, delimiter = '\t')
            for l in r:
                pe_ref[l['protein'].replace('NX_','')] = l['pe']

    protein_sequence = ""
    protein_name = ""
    gene_name = ""
    PE_number = "10"
    proteome = Proteome()
    raw_lines = []

    for line in open(filename, "r"):
        line = line.strip()
        if len(line) == 0:
            continue

        if line[0] == ">":
            if not len(protein_sequence) == 0:
                new_protein = Protein(protein_name, protein_sequence, gene_name)
                new_protein.pe_number = PE_number
                proteome.add_protein(new_protein)


            protein_name = line.split(" ")[0][1:].split(";")[0].replace("GFFtoFASTA@","")
            raw_lines.append(line)
            protein_sequence = ""
            #Try to get gene name
            name_splits = line.split(" ")
            for split in name_splits:
                if split.find("GN=") != -1:
                    gene_name = split[3:]
                if split.find("PE=") != -1:
                    PE_number = split[3:]
            if (external_PE):
                try:
                    PE_number = pe_ref.get(protein_name.split("|")[1],'0')
                except:
                    PE_number = '0'
        else:
            protein_sequence += line.rstrip()

    #Get last protein
    new_protein = Protein(protein_name, protein_sequence, gene_name)
    new_protein.pe_number = PE_number
    proteome.add_protein(new_protein)

    if output_raw:
        return proteome, raw_lines
    else:
        return proteome

class Proteome:
    def __init__(self):
        self.protein_list = []
        self.protein_map = {}
        self.fourmer_hash_to_proteins = None


    def write_summary(self, output_filename):
        ming_fileio_library.write_dictionary_table_data(self.produce_protein_dict(), output_filename)

    def produce_protein_dict(self):
        output_dict = defaultdict(list)

        for protein in self.protein_list:
            output_dict["protein"].append(protein.protein)
            output_dict["score"].append(protein.score)
            output_dict["number_precursors"].append(protein.number_precursors)
            output_dict["fdr"].append(protein.fdr)
            output_dict["pe_number"].append(protein.pe_number)
            output_dict["length"].append(len(protein.sequence))

        return output_dict


    def add_protein(self, new_protein):
        self.protein_list.append(new_protein)
        self.protein_map[new_protein.protein] = new_protein


    def get_protein(self, protein_name):
        if protein_name in self.protein_map:
            return self.protein_map[protein_name]
        return None
        #for protein in self.protein_list:
        #    if protein.protein == protein_name:
        #        return protein
        #return None

    def get_peptides_mapped_to_proteins(self, peptide_list):
        peptide_mapping = defaultdict(list)
        for peptide in peptide_list:
            for protein in self.protein_list:
                if protein.contain_peptide(peptide) == True:
                    peptide_mapping[peptide].append(protein.protein)
        return peptide_mapping







    #This is the fast way to do this
    def get_peptides_mapped_to_proteins_efficient(self, peptide_list):
        peptide_mapping = defaultdict(list)

        fourmer_hash_to_proteins = None
        if self.fourmer_hash_to_proteins == None:
            #hash 4-mers
            fourmer_hash_to_proteins = defaultdict(list)
            for protein in self.protein_list:
                substring_set = set()
                for i in range(len(protein.sequence)):
                    substring = protein.sequence[i:i+4]
                    substring_set.add(substring)
                for substring in substring_set:
                    fourmer_hash_to_proteins[substring].append(protein.protein)

            self.fourmer_hash_to_proteins = fourmer_hash_to_proteins
        else:
            fourmer_hash_to_proteins = self.fourmer_hash_to_proteins

        print("Finished Precomputing Mapping Hash")
        mapping_progress = 0

        for peptide in peptide_list:
            mapping_progress += 1
            protein_candidates_prefix = set(fourmer_hash_to_proteins[peptide[:4]])
            protein_candidates_suffix = set(fourmer_hash_to_proteins[peptide[-4:]])

            protein_candidates = protein_candidates_prefix & protein_candidates_suffix

            #print(len(protein_candidates_prefix), len(protein_candidates_suffix), len(protein_candidates))
            if mapping_progress % 1000 == 0:
                print(mapping_progress, "of", len(peptide_list))

            for protein_name in protein_candidates:
                if self.protein_map[protein_name].contain_peptide(peptide) == True:
                    peptide_mapping[peptide].append(protein_name)

        return peptide_mapping

    from collections import defaultdict, namedtuple

    def find_overlapping(self, peptide_len, prefixes, suffixes):
        overlapping = []
        prefixes = sorted(prefixes, key = lambda x: (x[1],x[0]))
        suffixes = sorted(suffixes, key = lambda x: (x[1],x[0]))
        prefix_pos = 0
        suffix_pos = 0
        while(prefix_pos < len(prefixes) and suffix_pos < len(suffixes)):
            if prefixes[prefix_pos][1] == suffixes[suffix_pos][1]:
                if prefixes[prefix_pos][0] + peptide_len - 4 == suffixes[suffix_pos][0]:
                    overlapping.append(prefixes[prefix_pos])
                    prefix_pos += 1
                    suffix_pos += 1
                elif prefixes[prefix_pos][0] + peptide_len - 4 < suffixes[suffix_pos][0]:
                    prefix_pos += 1
                else:
                    suffix_pos += 1
            else:
                if prefixes[prefix_pos][1] < suffixes[suffix_pos][1]:
                    prefix_pos += 1
                else:
                    suffix_pos += 1
        return overlapping


    def find_match(self, peptide, protein_map, prefixes, suffixes):
        overlapping = set()
        prefixes = [(True,prefix) for prefix in prefixes]
        suffixes = [(False,suffix) for suffix in suffixes]
        for is_prefix, (pos, protein) in prefixes + suffixes:
            if is_prefix:
                start_pos = pos
                end_pos = pos+len(peptide)
            else:
                start_pos = pos+4-len(peptide)
                end_pos = pos+4
            protein_subsequence = protein_map[protein].sequence[start_pos:end_pos]
            matches = 0
            il_matches = 0
            for aa_protein, aa_peptide in zip(protein_subsequence, peptide):
                if aa_protein.replace('I','L') == aa_peptide.replace('I','L'):
                    il_matches += 1
                    if aa_protein == aa_peptide:
                        matches += 1
            if il_matches >= len(peptide)-1:
                # should be one-based, don't need to adjust end because end is currently non-inclusive
                overlapping.add(ProteomeMatch(protein,start_pos+1,end_pos,matches!=il_matches,len(peptide)-il_matches))
        return overlapping

    #This is the fast way to do this
    def get_peptides_mapped_to_proteins_efficient_w_coordinates_w_mismatch(self, peptide_list):
        peptide_mapping = defaultdict(list)

        fourmer_hash_to_proteins = defaultdict(list)
        for protein in self.protein_list:
            substring_set = set()
            for i in range(len(protein.sequence)):
                #do the IL substitution before mapping, since we want to keep track of both exact matches and I/L matches
                substring = protein.sequence[i:i+4].replace('I','L')
                fourmer_hash_to_proteins[substring].append((i,protein.protein))

        for peptide in peptide_list:
            il_peptide = peptide.replace('I','L')

            protein_candidates_prefix = fourmer_hash_to_proteins[il_peptide[:4]]
            protein_candidates_suffix = fourmer_hash_to_proteins[il_peptide[-4:]]

            protein_candidates = self.find_match(peptide, self.protein_map, protein_candidates_prefix, protein_candidates_suffix)
            peptide_mapping[peptide] = list(protein_candidates)

        return peptide_mapping

    #This is the fast way to do this
    def get_peptides_mapped_to_proteins_efficient_w_coordinates(self, peptide_list):
        peptide_mapping = defaultdict(list)

        # fourmer_hash_to_proteins = None
        # if self.fourmer_hash_to_proteins == None:
            #hash 4-mers
        fourmer_hash_to_proteins = defaultdict(list)
        for protein in self.protein_list:
            substring_set = set()
            for i in range(len(protein.sequence)):
                substring = protein.sequence[i:i+4]
                fourmer_hash_to_proteins[substring].append((i,protein.protein))
                #substring_set.add(substring)
            #for substring in substring_set:

        self.fourmer_hash_to_proteins = fourmer_hash_to_proteins
        # else:
        #     fourmer_hash_to_proteins = self.fourmer_hash_to_proteins

        print("Finished Precomputing Mapping Hash")
        mapping_progress = 0

        for peptide in peptide_list:
            mapping_progress += 1
            protein_candidates_prefix = fourmer_hash_to_proteins[peptide[:4]]
            protein_candidates_suffix = fourmer_hash_to_proteins[peptide[-4:]]

            protein_candidates = self.find_overlapping(len(peptide),protein_candidates_prefix,protein_candidates_suffix)

            #print(len(protein_candidates_prefix), len(protein_candidates_suffix), len(protein_candidates))
            if mapping_progress % 1000 == 0:
                print(mapping_progress, "of", len(peptide_list))
            if len(protein_candidates) == 0:
                continue

            for p in protein_candidates:
                pos = p[0]
                protein_name = p[1]
                if self.protein_map[protein_name].sequence[pos:pos+len(peptide)] == peptide:
                    peptide_mapping[peptide].append("{} ({}-{})".format(protein_name,pos,pos+len(peptide)))
        return peptide_mapping



    #Does not match short peptides less than length 9
    def get_peptides_mapped_to_proteins_mismatched(self, peptide_list, mismatch_max=1):
        peptide_mapping = defaultdict(list)

        fourmer_hash_to_proteins = None
        if self.fourmer_hash_to_proteins == None:
            #hash 4-mers
            fourmer_hash_to_proteins = defaultdict(list)
            for protein in self.protein_list:
                substring_set = set()
                for i in range(len(protein.sequence)):
                    substring = protein.sequence[i:i+4]
                    substring_set.add(substring)
                for substring in substring_set:
                    fourmer_hash_to_proteins[substring].append(protein.protein)

            self.fourmer_hash_to_proteins = fourmer_hash_to_proteins
        else:
            fourmer_hash_to_proteins = self.fourmer_hash_to_proteins

        print("Finished Precomputing Mapping Hash")
        mapping_progress = 0

        for peptide in peptide_list:
            mapping_progress += 1

            protein_candidates = set()

            if len(peptide) < 8:
                continue

            if len(peptide) < 12:
                protein_candidates_prefix = set(fourmer_hash_to_proteins[peptide[:4]])
                protein_candidates_suffix = set(fourmer_hash_to_proteins[peptide[-4:]])

                protein_candidates = protein_candidates_prefix.union(protein_candidates_suffix)

            if len(peptide) >= 12:
                protein_candidates_prefix = set(fourmer_hash_to_proteins[peptide[:4]])
                protein_candidates_suffix = set(fourmer_hash_to_proteins[peptide[-4:]])
                protein_candidates_middle = set(fourmer_hash_to_proteins[peptide[4:8]])

                intersection1 = protein_candidates_prefix.intersection(protein_candidates_suffix)
                intersection2 = protein_candidates_prefix.intersection(protein_candidates_middle)
                intersection3 = protein_candidates_suffix.intersection(protein_candidates_middle)

                protein_candidates = intersection1 | intersection2 | intersection3

            #print(len(protein_candidates_prefix), len(protein_candidates_suffix), len(protein_candidates))
            if mapping_progress % 1000 == 0:
                print(mapping_progress, "of", len(peptide_list))

            #print(mapping_progress, "of", len(peptide_list), len(protein_candidates), peptide)

            for protein_name in protein_candidates:
                if self.protein_map[protein_name].contain_peptide_mismatch(peptide, mismatch_max=1) == True:
                    peptide_mapping[peptide].append(protein_name)

        return peptide_mapping

    #
    # def get_peptides_mapped_to_proteins_bwt(self, peptide_list, mismatch_max=0):
    #     peptide_mapping = defaultdict(list)
    #
    #     fourmer_hash_to_proteins = None
    #     if self.fourmer_hash_to_proteins == None:
    #         #hash 4-mers
    #         fourmer_hash_to_proteins = defaultdict(list)
    #         for protein in self.protein_list:
    #             substring_set = set()
    #             for i in range(len(protein.sequence)):
    #                 substring = protein.sequence[i:i+4]
    #                 substring_set.add(substring)
    #             for substring in substring_set:
    #                 fourmer_hash_to_proteins[substring].append(protein.protein)
    #
    #         self.fourmer_hash_to_proteins = fourmer_hash_to_proteins
    #     else:
    #         fourmer_hash_to_proteins = self.fourmer_hash_to_proteins
    #
    #     print("Finished Precomputing Mapping Hash")
    #     protein_count = 0
    #     for protein in self.protein_list:
    #         protein_count += 1
    #         #Debugging
    #         #if protein_count > 1000:
    #         #    continue
    #         #print(protein.protein, protein_count)
    #         bwt_data = bwt.make_all(protein.sequence)
    #         protein.bwt_data = bwt_data
    #
    #     print("Finished Precomputing BWT")
    #     mapping_progress = 0
    #     for peptide in peptide_list:
    #         mapping_progress += 1
    #         protein_candidates_prefix = set(fourmer_hash_to_proteins[peptide[:4]])
    #         protein_candidates_suffix = set(fourmer_hash_to_proteins[peptide[-4:]])
    #
    #         protein_candidates = protein_candidates_prefix.union(protein_candidates_suffix)
    #
    #         if mapping_progress % 1000 == 0:
    #             print(mapping_progress, "of", len(peptide_list))
    #
    #         for protein_name in protein_candidates:
    #             if "bwt_data" in self.protein_map[protein_name].__dict__:
    #                 found_locations = bwt.find(peptide, self.protein_map[protein_name].sequence, mismatches=mismatch_max, bwt_data=self.protein_map[protein_name].bwt_data)
    #                 if len(found_locations) > 0:
    #                     #print(peptide, protein_name, found_locations)
    #                     peptide_mapping[peptide].append(protein_name)
    #
    #     return peptide_mapping






    def get_proteins_to_peptides_efficient(self, peptide_list):
        protein_mapping = defaultdict(list)

        #hash 4-mers
        fourmer_hash_to_proteins = defaultdict(list)
        for protein in self.protein_list:
            substring_set = set()
            for i in range(len(protein.sequence)):
                substring = protein.sequence[i:i+4]
                substring_set.add(substring)
            for substring in substring_set:
                fourmer_hash_to_proteins[substring].append(protein.protein)

        for peptide in peptide_list:
            protein_candidates = fourmer_hash_to_proteins[peptide[:4]]
            for protein_name in protein_candidates:
                if self.protein_map[protein_name].contain_peptide(peptide) == True:
                    protein_mapping[protein_name].append(peptide)
        return protein_mapping

    #Return a list of proteins that are covered by a set of peptides
    #Not really teh best way to do this
    #DEPRECATED
    def get_proteins_covered_by_peptides_hash(self, peptide_list, do_reverse=False, protein_pkl_filename = None):
        peptide_to_protein_hash = {}
        create_hash = False

        if protein_pkl_filename == None:
            create_hash = True
        else:
            if os.path.exists(protein_pkl_filename):
                #Load that shit
                create_hash = False
                peptide_to_protein_hash = pickle.load(open(protein_pkl_filename, 'rb'))
            else:
                create_hash = True


        if create_hash == True:
            print("Creating Hash")
            protein_count = 0
            for protein in self.protein_list:
                protein_count += 1
                print(protein.protein + "\t" + str(protein_count) + " of " + str(len(self.protein_list)))
                #print peptide_to_protein_hash
                for i in range(len(protein.sequence)):
                    for j in range(len(protein.sequence)):
                        if j < i or j - i > 20 or j - i < 4:
                            continue
                        substring = protein.sequence[i:j+1]
                        if not substring in peptide_to_protein_hash:
                            peptide_to_protein_hash[substring] = set()
                        peptide_to_protein_hash[substring].add(protein.protein)

            if protein_pkl_filename != None:
                pickle.dump(peptide_to_protein_hash, open(protein_pkl_filename, 'rb'))

        #Doing Actual Search
        search_peptide_list = []

        if do_reverse == True:
            for peptide in peptide_list:
                rev_peptide = peptide[::-1]
                search_peptide_list.append(peptide)
                search_peptide_list.append(rev_peptide)
        else:
            search_peptide_list = peptide_list

        #Do shit now with it
        return_proteins = []
        for peptide in search_peptide_list:
            if peptide in peptide_to_protein_hash:
                return_proteins += list(peptide_to_protein_hash[peptide])

        return return_proteins


    #Return a list of proteins that are covered by a set of peptides
    def get_proteins_covered_by_peptides(self, peptide_list, do_reverse=False):
        return_proteins = []

        search_peptide_list = []

        if do_reverse == True:
            for peptide in peptide_list:
                rev_peptide = peptide[::-1]
                search_peptide_list.append(peptide)
                search_peptide_list.append(rev_peptide)
        else:
            search_peptide_list = peptide_list

        #Lets make one big string and see how it goes
        protein_sequence_list = []
        protein_position_list = []
        protein_name_list = []
        running_position = 0
        for protein in self.protein_list:
            protein_sequence_list.append(protein.sequence)
            protein_position_list.append(running_position)
            protein_name_list.append(protein)
            running_position += len(protein.sequence) + 1

        #make big string
        proteome_string = ";".join(protein_sequence_list)

        peptide_count = 0
        for peptide in search_peptide_list:
            peptide_count += 1
            print(str(peptide_count) + " of " + str(len(search_peptide_list)))

            find_position = 0
            while find_position != -1:
                find_position = proteome_string.find(peptide, find_position + 1)

                if find_position == -1:
                    break

                #Determine which protein is that position
                for i in range(len(protein_position_list)):
                    if find_position < protein_position_list[i]:
                        protein_found = protein_name_list[i-1]
                        return_proteins.append(protein_found)
                        break

        return return_proteins

    #Return a list of proteins that are covered by a set of peptides with at least K per protein
    def get_proteins_covered_by_k_peptides(self, peptide_list, minimum_peptides=1, do_reverse=False):
        return_proteins = []
        search_peptide_list = []

        if do_reverse == True:
            for peptide in peptide_list:
                rev_peptide = peptide[::-1]
                search_peptide_list.append(peptide)
                search_peptide_list.append(rev_peptide)
        else:
            search_peptide_list = peptide_list

        search_peptide_list = list(set(search_peptide_list))

        protein_name_to_number_of_peptides_count = self.get_proteins_with_number_of_peptides_covered_map(search_peptide_list)

        for protein in self.protein_list:
            number_of_peptides = protein_name_to_number_of_peptides_count[protein.protein]
            if number_of_peptides >= minimum_peptides:
                return_proteins.append(protein)

        return return_proteins

    #Given a list of peptides, returns a map of protein names to the number of peptides that cover it
    def get_proteins_with_number_of_peptides_covered_map(self, peptide_list):
        protein_to_peptides_map = self.get_proteins_to_peptides_efficient(peptide_list)

        return_proteins = {}
        for protein in self.protein_list:
            number_of_peptides = 0

            if protein.protein in protein_to_peptides_map:
                number_of_peptides = len(protein_to_peptides_map[protein.protein])

            return_proteins[protein.protein] = number_of_peptides
        return return_proteins

    #List of peptide sequences
    def calculate_peptide_coverage_of_all_proteins(self, peptides):
        coverage_map = {}
        count = 0
        for protein in self.protein_list:
            count += 1
            coverage_number = protein.coverage(peptides)
            coverage_map[protein.protein] = coverage_number

            #print str(count) + " of " + str(len(self.protein_list))
            #if coverage_number > 0:
            #    print protein.protein + "\t" + str(coverage_number)
        return coverage_map

    #Prints out the percetnage of tryptic peptides covered per protein
    def calculate_tryptic_peptide_coverage_of_all_proteins(self, peptides):
        coverage_map = {}
        for protein in self.protein_list:
            coverage_number = protein.get_tryptic_peptide_coverage(peptides)
            coverage_map[protein.protein] = coverage_number
            #if coverage_number > 0:
            #    print protein.protein + "\t" + str(coverage_number)
        return coverage_map

    def calculate_unique_trypic_peptide_coverage_of_all_proteins(self, peptides):
        coverage_map = {}
        count = 0
        unique_peptides = self.get_unique_peptides()
        for protein in self.protein_list:
            count += 1
            #print str(count) + " of " + str(len(self.protein_list))
            coverage_number = protein.get_unique_tryptic_peptide_coverage(peptides, unique_peptides)
            coverage_map[protein.protein] = coverage_number
            #if coverage_number > 0:
            #    print protein.protein + "\t" + str(coverage_number)
        return coverage_map

    #Returns a list of protein objects for where the peptide matches
    def get_proteins_with_sequence(self, peptide):
        protein_list = []
        for protein in self.protein_list:
            if protein.contain_peptide(peptide):
                protein_list.append(protein)
        return protein_list

    #Returns a list of all trypic peptides
    def get_tryptic_peptides(self, min_length=0):
        tryptic_peptides = []
        for protein in self.protein_list:
            tryptic_peptides += protein.get_tryptic_peptides(min_length)

        return list(set(tryptic_peptides))

    #Given a list of peptides that were identified
    #Returns the percentage covered of proteome, percentage of input peptides in proteome
    def calculate_tryptic_peptides_covered(self, covered_peptides):
        proteome_peptides = self.get_tryptic_peptides()
        #proteome_peptides_map = {}
        #output_proteome_file = open("output_proteome.out", "w")
        #for peptide in proteome_peptides:
        #    proteome_peptides_map[peptide] = 1
        #    output_proteome_file.write(peptide + "\n")

        #intersection_set = []
        #for peptide in covered_peptides:
        #    if peptide in proteome_peptides_map:
        #        intersection_set.append(peptide)

        intersection_set = (set(covered_peptides)).intersection(set(proteome_peptides))
        return float(len(intersection_set))/float(len(proteome_peptides)), float(len(intersection_set))/float(len(covered_peptides))

    def set_coverage_with_proteome(self, covered_peptides):
        proteome_peptides = self.get_tryptic_peptides()
        intersection_set = list((set(covered_peptides)).intersection(set(proteome_peptides)))
        cover_exclusive = list((set(covered_peptides)).difference(set(proteome_peptides)))
        proteome_exclusive = list((set(proteome_peptides)).difference(set(covered_peptides)))
        return intersection_set, proteome_exclusive, cover_exclusive

    #Give set of unique peptides that are unique to proteins
    def get_unique_peptides(self, min_length=0):
        peptide_sharing_count = {}
        pattern = re.compile('[KR][^P]')
        for protein in self.protein_list:
            tryptic_peptides = list(set(protein.get_tryptic_peptides(min_length)))
            for peptide in tryptic_peptides:
                if not peptide in peptide_sharing_count:
                    peptide_sharing_count[peptide] = 0
                peptide_sharing_count[peptide]+= 1

        unique_peptides = []
        for peptide in peptide_sharing_count:
            if peptide_sharing_count[peptide] == 1:
                unique_peptides.append(peptide)

        return unique_peptides

    def calculate_unique_peptide_coverage(self, covered_peptides):
        unique_peptides = self.get_unique_peptides()
        intersection_set = (set(covered_peptides)).intersection(set(unique_peptides))
        return float(len(intersection_set))/float(len(unique_peptides))

    #Get Genes Covered
    def get_genes_covered_by_peptides(self, peptides):
        genes_covered = []

        for protein in self.protein_list:
            coverage_number = protein.coverage(peptides)
            if coverage_number > 0:
                genes_covered.append(protein.gene_name)
        return list(set(genes_covered))

    #Get All Genes
    def get_all_genes(self):
        all_genes = []
        for protein in self.protein_list:
            all_genes.append(protein.gene_name)
        return list(set(all_genes))

    #Returns a list of peptides that uniquely identify a gene
    def get_unique_gene_peptides(self):
        gene_to_peptide = {} #Map of all peptides per gene
        for protein in self.protein_list:
            if not protein.gene_name in gene_to_peptide:
                gene_to_peptide[protein.gene_name] = []
            gene_to_peptide[protein.gene_name] += protein.get_tryptic_peptides()

        #Making the peptides unique
        for gene in gene_to_peptide:
            gene_to_peptide[gene] = list(set(gene_to_peptide[gene]))

        #Count Uniqueness of peptides
        peptide_sharing_count = {}
        for gene in gene_to_peptide:
            for peptide in gene_to_peptide[gene]:
                if not peptide in peptide_sharing_count:
                    peptide_sharing_count[peptide] = 0
                peptide_sharing_count[peptide] += 1

        #Find Unique Peptides
        unique_peptides = []
        for peptide in peptide_sharing_count:
            if peptide_sharing_count[peptide] == 1:
                unique_peptides.append(peptide)

        return unique_peptides

    def create_decoy_proteins(self):
        decoy_proteins = []
        for protein in self.protein_list:
            decoy_protein = Protein("XXX_" + protein.protein, protein.sequence[::-1], protein.gene_name)
            decoy_protein.decoy = 1
            decoy_proteins.append(decoy_protein)

        for protein in decoy_proteins:
            self.add_protein(protein)

    def get_proteome_amino_acid_count(self):
        total_amino_acids = 0
        for protein in self.protein_list:
            total_amino_acids += len(protein.sequence)
        return total_amino_acids


class Protein:
    def __init__(self, protein_name, sequence, gene_name):
        self.protein = protein_name
        self.sequence = sequence
        self.peptides = []
        self.gene_name = gene_name
        self.decoy = 0
        self.pe_number = ""
        self.fdr = -1.0
        self.precursors_mapped = 0
        self.score = 0


    #Returns the coverage of protein, given peptides
    def coverage(self, peptides):
        coverage_list = self.coverage_by_amino_acids(peptides)
        return float(sum(coverage_list))/float(len(coverage_list))

    def coverage_by_amino_acids(self, peptides):
        coverage_list = [0] * len(self.sequence)

        if len(coverage_list) == 0:
            return 0

        for peptide in peptides:
            starting_location = 0
            found_location = self.sequence.find(peptide, starting_location)
            #print self.sequence
            while found_location != -1:
                #coverage_list[found_location:found_location + len(self.sequence)] = 1
                for i in range(found_location, found_location + len(peptide)):
                    coverage_list[i] = 1

                starting_location = found_location + 1
                found_location = self.sequence.find(peptide, starting_location)

        return coverage_list

    def coverage_of_uncoverable_sequences(self):
        tryptic_peptides = self.get_tryptic_peptides(min_length=4)
        unobservable_peptides = []
        for peptide in tryptic_peptides:
            #if len(peptide) < 6 or len(peptide) > 40:
            if len(peptide) > 40:
                unobservable_peptides.append(peptide)

        coverage_list = self.coverage_by_amino_acids(unobservable_peptides)

        return coverage_list

    def coverage_of_coverable_sequences(self):
        uncoverable_list = self.coverage_of_uncoverable_sequences()
        coverable_list = []
        for coverage_site in uncoverable_list:
            if coverage_site == 0:
                coverable_list.append(1)
            else:
                coverable_list.append(0)

        return coverable_list






    def number_of_non_overlapping_sequences(self, peptides):
        number_of_non_overlapping_sequences = 0

        coverage_list = [0] * len(self.sequence)
        if len(coverage_list) == 0:
            return 0

        for peptide in peptides:
            starting_location = 0
            found_location = self.sequence.find(peptide, starting_location)

            all_locations = []

            while found_location != -1:
                #coverage_list[found_location:found_location + len(self.sequence)] = 1
                for i in range(found_location, found_location + len(peptide)):
                    all_locations.append(i)

                starting_location = found_location + 1
                found_location = self.sequence.find(peptide, starting_location)

            #check if we should add this peptide
            not_overlapping = True
            for location in all_locations:
                if coverage_list[location] == 1:
                    not_overlapping = False
                    break

            if not_overlapping == True:
                for location in all_locations:
                    coverage_list[location] = 1
                number_of_non_overlapping_sequences += 1

        return number_of_non_overlapping_sequences





    def map_location_of_peptide(self, peptide):
        found_location = self.sequence.find(peptide)
        return found_location

    #returns percentage of tryptic peptides covered in protein
    def get_tryptic_peptide_coverage(self, peptides):
        tryptic_peptides = set(self.get_tryptic_peptides())
        found_tryptic_peptides = set()

        for peptide in peptides:
            if peptide in tryptic_peptides:
                found_tryptic_peptides.add(peptide)

        return float(len(found_tryptic_peptides))/float(len(tryptic_peptides))

    def get_unique_tryptic_peptide_coverage(self, peptides, unique_peptides):
        tryptic_peptides = set(self.get_tryptic_peptides())
        unique_tryptic_peptides = set(unique_peptides).intersection(tryptic_peptides)
        found_tryptic_peptides = set()

        if(len(unique_tryptic_peptides) == 0):
            return -1

        for peptide in peptides:
            if peptide in unique_tryptic_peptides:
                found_tryptic_peptides.add(peptide)

        return float(len(found_tryptic_peptides))/float(len(unique_tryptic_peptides))

    def contain_peptide(self, peptide):
        if peptide in self.sequence:
            return True
        else:
            return False

    def contain_peptide_mismatch(self, peptide, mismatch_max=1):
        #length_of_query = len(peptide)
        #for i in range(len(self.sequence) - length_of_query):
        #    excised_peptide = self.sequence[i:i+length_of_query]
        #    hamming_distance = distance.hamming(peptide, excised_peptide)
        #    if hamming_distance <= mismatch_max:
        #        return True
        #return False
        regex_string = "(%s){e<=1}" % (peptide)
        m=regex.findall(regex_string, self.sequence)
        if len(m) > 0:
            return True
        return False



    def number_of_peptides_on_protein(self, peptides):
        peptide_map_count = 0
        for peptide in peptides:
            if self.contain_peptide(peptide):
                peptide_map_count += 1
        return peptide_map_count

    #Tryptic Peptides for protein
    def get_tryptic_peptides(self, min_length=0):
        tryptic_peptides = re.findall(r".(?:(?<![KR](?!P)).)*", self.sequence)
        return_list = []
        for peptide in tryptic_peptides:
            if len(peptide) > min_length:
                return_list.append(peptide)
        return return_list

    def to_fasta(self):
        output_string = ">" + self.protein +  " " + "GN=" + self.gene_name + "\n"
        output_string += self.sequence.strip()
        return output_string
