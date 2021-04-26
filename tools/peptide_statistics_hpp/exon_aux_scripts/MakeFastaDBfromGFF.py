'''
Created on Jun 16, 2017

@author: s3cha
'''

import os
import pysam
import sys
import ming_proteosafe_library
import json


class GTFtoFASTA(object):
    write_stopcodon = 0
    DNA = {}
    fasta_index = 0
    ForwardCode = {"TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
               "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",
               "TAT":"Y", "TAC":"Y", "TAA":"X", "TAG":"X",
               "TGT":"C", "TGC":"C", "TGA":"X", "TGG":"W",
               "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
               "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
               "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
               "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
               "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
               "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
               "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
               "AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",
               "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
               "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
               "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
               "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G",
               }
    ReverseCode =   {"AAA":"F", "AAG":"F", "AAT":"L", "AAC":"L",
               "AGA":"S", "AGG":"S", "AGT":"S", "AGC":"S",
               "ATA":"Y", "ATG":"Y", "ATT":"X", "ATC":"X",
               "ACA":"C", "ACG":"C", "ACT":"X", "ACC":"W",
               "GAA":"L", "GAG":"L", "GAT":"L", "GAC":"L",
               "GGA":"P", "GGG":"P", "GGT":"P", "GGC":"P",
               "GTA":"H", "GTG":"H", "GTT":"Q", "GTC":"Q",
               "GCA":"R", "GCG":"R", "GCT":"R", "GCC":"R",
               "TAA":"I", "TAG":"I", "TAT":"I", "TAC":"M",
               "TGA":"T", "TGG":"T", "TGT":"T", "TGC":"T",
               "TTA":"N", "TTG":"N", "TTT":"K", "TTC":"K",
               "TCA":"S", "TCG":"S", "TCT":"R", "TCC":"R",
               "CAA":"V", "CAG":"V", "CAT":"V", "CAC":"V",
               "CGA":"A", "CGG":"A", "CGT":"A", "CGC":"A",
               "CTA":"D", "CTG":"D", "CTT":"E", "CTC":"E",
               "CCA":"G", "CCG":"G", "CCT":"G", "CCC":"G",
               }

    def __init__(self):
        pass

    def read_dna_trie(self,folder,params):
        files = os.listdir(folder)
        for file in files:
            if not file.endswith('.trie'):
                continue
            chromosome = params[file].split('/')[-1].split('_')[0]
            print(chromosome)
            self.DNA[chromosome] = open(folder.rstrip('/')+'/'+file,'r').readline()
#             break
        return self.DNA

    def hamming_distance(self,s1,s2):
        if len(s1) != len(s2):
            return sys.maxint
        distance = 0
        for index in range(len(s1)):
            if s1[index] != s2[index]:
                distance += 1
        return distance

    def getAA(self,output_array,strand):
        if strand == '-':
            output_array = output_array[::-1]
        AA_seq = ''
        for i in range(len(output_array[:])/3):
            if 'N' in output_array[3*i:3*(i+1)]:
                AA_seq += 'X'
            else:
                if strand == '-':
                    AA_seq += self.ReverseCode[output_array[3*i:3*(i+1)]]
                else:
                    AA_seq += self.ForwardCode[output_array[3*i:3*(i+1)]]
        return AA_seq

    def convertBase_one2zero(self,exons):
        return [(exon[0]-1,exon[1]) for exon in exons]

    def WriteFASTA(self,chr,exons,strand,s,genename=None,geneid=None,frame_error = 0):
        if not exons:
            return
        if chr not in self.DNA:
            return
        if frame_error == 0:
            exons = self.convertBase_one2zero(exons)
            if strand == '+':
                exons[-1] = (exons[-1][0],exons[-1][1]-3)
            else:
                exons[-1] = (exons[-1][0]+3,exons[-1][1])

        def write_header(current_exons,temp_strand):
            temp_strand = 0 if temp_strand == '-' else 1
            exon_coord = ';'.join(['/'.join(map(str,exon)) for exon in current_exons])+';'
            s.write('>GFFtoFASTA@%s@%s@%s@%s;%s\n'%(chr,genename,geneid,temp_strand,exon_coord))
            self.fasta_index += 1
            return

        def write_sequence(seq,temp_strand):
            s.write('%s\n'%(self.getAA(seq,temp_strand)))
            return

        for index_frame in range(1):
            index_frame = frame_error
            length = sum([exon[1]-exon[0] for exon in exons]) - index_frame
            trim_length = length%3
            if strand == '+':

                current_exons = [(exons[0][0]+index_frame,exons[0][1])]
                current_exons += [exon for exon in exons[1:len(exons)-1]]
                if len(exons) > 1:
                    current_exons.append((exons[-1][0],exons[-1][1]-trim_length))
                seq = ''.join([self.DNA[chr][exon[0]:exon[1]] for exon in current_exons])
            else:

                current_exons = [(exons[0][0],exons[0][1]-index_frame)]
                current_exons += [exon for exon in exons[1:len(exons)-1]]
                if len(exons) > 1:
                    current_exons.append((exons[-1][0]+trim_length,exons[-1][1]))
                seq = ''.join([self.DNA[chr][exon[0]:exon[1]] for exon in current_exons[::-1]])
            if len(seq.strip()) <= 9:
                continue
            if strand != '-':
                aaseq = self.getAA(seq,'+')
                if 'X' in aaseq:
                    if frame_error > 2:
                        print 'Frame error: ',strand,genename,geneid, frame_error,aaseq
                        write_header(current_exons,'+')
                        write_sequence(seq,'+')
                        return
                    self.WriteFASTA(chr,exons,strand,s,genename,geneid,frame_error +1)
                    continue
                write_header(current_exons,'+')
                write_sequence(seq,'+')
            if strand != '+':
                aaseq = self.getAA(seq,'-')

                if 'X' in aaseq:
                    if frame_error > 2:
                        print 'Frame error: ',strand,genename,geneid, frame_error,aaseq
                        write_header(current_exons,'-')
                        write_sequence(seq,'-')
                        return
                    self.WriteFASTA(chr,exons,strand,s,genename,geneid,frame_error +1)
                    continue
                write_header(current_exons,'-')
                write_sequence(seq,'-')
#         print chr,exons,strand
#         if geneid == 'Parent=ENSG00000165322,Name=ARHGAP12-005':
#             sys.exit()
#         print aaseq
#         sys.exit()

        return

    def GFFtoFASTA(self,gff,output,transcript_identifier = 'mRNA'):
        s = open(output,'w')
        if not self.DNA:
            raise ValueError('No DNA file has been included')
        gff_header = []
        exons = []

        count = 0
        transcriptid = None
        geneid = ''
        genename = ''
        chr = self.DNA.keys()[0]
        strand = ''
        for line in open(gff,'r'):
#             if count > 1000:
#                 break
            if line.startswith('#'):
                gff_header.append(line.strip())
                continue
            line = line.strip().split('\t')
            if len(line) < 7:
                continue
            chr, type, start, end = [line[0],line[2],int(line[3]),int(line[4])]
            if chr not in self.DNA:
                continue
#             if type == 'gene':
#                 info = line[8].strip().split(';')
#                 genename = info[1]
#                 geneid = info[0]
            if type == transcript_identifier:
                if exons:
                    self.WriteFASTA(chr,exons,strand,s,"EN|" + transcriptid,"EN|" + geneid+','+genename)
                info = line[8].strip().split(';')
                if '"' not in info[0]:
                    transcriptid = info[0]
                    geneid =  info[2]
                    genename = info[1]
                    chr = line[0]
                    strand = line[6]
                    exons = []
                    continue
                else:
                    temp = {}
                    for part_info in info:
                        part_info = part_info.split('"')
                        if len(part_info)>1:
                            temp[part_info[0].strip()] = part_info[1]
                    transcriptid = temp.get('transcript_id')
                    geneid = temp.get('gene_id')
                    genename = temp.get('gene_name')
                    chr = line[0]
                    strand = line[6]
                    exons = []
                    continue
            elif type == 'CDS':
                exons.append((start,end))
            else:
#                 print 'Type %s is not defined'%(type)
                continue

            count += 1
        self.WriteFASTA(chr,exons,strand,s,"EN|" + transcriptid,"EN|" + geneid+','+genename)
#         self.WriteFASTA(chr,exons,strand,s,transcriptid,geneid)

    def GFFtoFASTA_temp(self,gff,output):
        s = open(output,'w')
        if not self.DNA:
            raise ValueError('No DNA file has been included')
        gff_header = []
        exons = []

        count = 0
        transcriptid = None
        geneid = ''
        genename = ''
        for line in open(gff,'r'):
            if line.startswith('#'):
                gff_header.append(line.strip())
                continue
            line = line.strip().split('\t')
            if len(line) < 7:
                continue
            chr, type, start, end, strand = [line[0],line[2],int(line[3]),int(line[4]),line[6]]
            if chr not in self.DNA:
                continue
#             if type == 'gene':
#                 info = line[8].strip().split(';')
#                 genename = info[1]
#                 geneid = info[0]
            if type == 'mRNA':
                count += 1
#                 print count,chr,exons
                self.WriteFASTA(chr,exons,strand,s,transcriptid,geneid+','+genename)
                info = line[8].strip().split(';')
                transcriptid = info[0]
                geneid =  info[2]
                genename = info[1]
                exons = []
                continue
            elif type == 'exon':
                exons.append((start,end))
            else:
#                 print 'Type %s is not defined'%(type)
                continue

        print count
        self.WriteFASTA(chr,exons,strand,s,transcriptid,geneid+','+genename)

    def remove_known_seq(self,knowndb,noveldb,newdb):
        f_known = open(knowndb,'r')
        known_db_seq = ''.join([line.strip() for line in f_known.readlines()])
        f_novel = open(noveldb,'r')
        s = open(newdb,'w')
        for line in f_novel:
            if line.startswith('>'):
                header = line
                continue
            else:
                if line.strip() not in known_db_seq:
                    s.write(header)
                    s.write(line)




if __name__ == '__main__':
    if len(sys.argv) < 4:
        version = 'GRCh38'
#         version = 'Uniprot'

        if version == 'GRCh37':
            dna_folder = '/home/s3cha/data/dna/GRCh37'
            gff_file = '/home/s3cha/data/Nuno/MakeAllTrypticPeptides/GFF/GRCh37.70/Homo_sapiens.GRCh37.70_formatted.gff'
            output_file = '/home/s3cha/data/Nuno/MakeAllTrypticPeptides/AlltrypticInEnsemble_GRCh37_CDS_data.fa'
        elif version == 'GRCh38':
            dna_folder = '/home/s3cha/data/dna/GRCh38'
    #         gff_file = '/home/s3cha/data/Nuno/MakeAllTrypticPeptides/GFF/GRCh38.p5/interim_GRCh38.p5_top_level_2015-09-25.gff3'
            gff_file = '/home/s3cha/data/Nuno/MakeAllTrypticPeptides/GFF/GRCh38.p5/Homo_sapiens.GRCh38.89.chr_reformatted.gtf'
            output_file = '/home/s3cha/data/Nuno/MakeAllTrypticPeptides/GRCh38_data/AlltrypticInEnsemble_GRCh38_CDS.fa'
        elif version == 'Uniprot':
            dna_folder = '/home/s3cha/data/dna/GRCh38'
    #         gff_file = '/home/s3cha/data/Nuno/MakeAllTrypticPeptides/GFF/GRCh38.p5/interim_GRCh38.p5_top_level_2015-09-25.gff3'
            gff_file = '/home/s3cha/data/Nuno/MakeAllTrypticPeptides/Uniprot_GRCh38_UCSCdown/Uniprot_GRCh38_UCSC_download.gtf'
            output_file = '/home/s3cha/data/Nuno/MakeAllTrypticPeptides/Uniprot_GRCh38_UCSCdown/Uniprot_GRCh38_CDS.fa'
    else:
#         print 'Input param is not initiated'
#         sys.exit()
        dna_folder = sys.argv[1]
        gff_file = sys.argv[2]
        output_file = sys.argv[3]
        input_parameters_filename = sys.argv[4]

    x = GTFtoFASTA()
    with open(input_parameters_filename) as f:
        params_map = ming_proteosafe_library.get_mangled_file_mapping(ming_proteosafe_library.parse_xml_file(f))
    dna = x.read_dna_trie(dna_folder, params_map)
    print dna.keys(),len(dna[dna.keys()[0]])
    x.GFFtoFASTA(gff_file,output_file,'transcript')
#     x.GFFtoFASTA(gff_file,output_file,'start_codon')

#     x.remove_known_seq(known_file,output_file,remove_known_output_file)
