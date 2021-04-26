'''
Created on Feb 28, 2017

@author: s3cha
'''
'''
Created on 2013. 3. 19.

@author: Akard3
'''
import os
import re
import sys
import time
from TagSearch import TagSearch
from collections import defaultdict
from csv import DictReader

def map_peptides_to_exons(exon_fasta, input_sequences):
    pepdic = {}
    for sequence in input_sequences:
        pepdic[sequence] = [[],{},0,0,0]

    ##data base read######################################################
    print('Reading database: ',exon_fasta)
    with open(exon_fasta,'r') as fasta_database:
        # fasta = fasta_database.readlines()
        fasta_info = []
        index_info = [0]
        fasta_seq = ''
        total_len = 0
        total_seq = []
        for count,fasta_line in enumerate(fasta_database):
            if count % 2 == 0:
                fasta_info.append(fasta_line)
            else:
                #sequence_info.append(fasta[i].strip())
                total_seq.append(fasta_line.strip().replace('I','L')+'X')
                total_len += len(fasta_line.strip().replace('I','L')+'X')
                index_info.append(total_len)
        fasta_seq = ''.join(total_seq)

        #>Splice@chr1@176@15444@0;20082126-20082213;20073655-20073753;20072948-20073092;20072024-20072144;20070131-20070219;  example of splice_info

        pep_list = pepdic.keys()
        count = 0
        tagS = TagSearch()
        tagS.makeTagTable(fasta_seq,4)
        for pep in pep_list:
            count += 1
            if count %1000 == 0:
                print count,':',time.ctime()
            original_pep = pep
            pep = CleanPeptideString(pep).replace('I','L')
    #         location_index = [m.start() for m in re.finditer(pep,fasta_seq)]
    #         location_index = [fasta_seq.find(pep)]
            location_index = tagS.getAllLocation(pep,fasta_seq)
            pep_location = []
            for location in location_index:
                index = binary_search(index_info,location,0,len(index_info)-1)
                splice_info = fasta_info[index].split('@')
                #full_seq = sequence_info[index]
                #start_seq = full_seq.find(pep)
                start_seq = location - index_info[index]
                end_seq = start_seq + len(pep)
                start_seq = start_seq *3
                end_seq = end_seq *3
                transcript = None
                gene = None
                if splice_info[0].find('Splice')>-1 or splice_info[0].find('Varient')>-1 or splice_info[0].find('rs')>-1 or splice_info[0].find('GFFtoFASTA') > -1:
                    chrNum = splice_info[1]
                    transcript = splice_info[2]
                    gene = splice_info[3]
                    splice_in = splice_info[4].split(';')
                    #print start_seq,end_seq, sequence, full_seq,
                    start = []
                    end = []
                    length = []
                    for i,splice in enumerate(splice_in):
                        if i == 0:
                            strand = int(splice)
                        elif i == len(splice_in)-1:
                            continue
                        else:
                            if splice.find('/')>-1:
                                splice = splice.split('/')
                            else:
                                splice = splice.split('-')
                            start.append(int(splice[0]))
                            end.append(int(splice[1]))
                            length.append(int(splice[1])-int(splice[0]))
                else:
                    if len(splice_info)<4:
                        continue
                    strand = int(splice_info[3])
                    start = [int(splice_info[1])]
                    end = [int(splice_info[2])]
                    length = [end[0]-start[0]]
                    chrNum = splice_info[0][1:]

                pep_location.append(get_location(start_seq,end_seq,start,end,length,strand,chrNum,gene,transcript))

            # prev = 0
            # i = 0
            # while i < len(pep_location):  # delete same location seq
            #     if prev == pep_location[i][3]:
            #         del pep_location[i]
            #     else:
            #         prev = pep_location[i][3]
            #         i += 1

            ### save pep location in pepdic with checking the occurence in previous pepdic
            for loc in pep_location:
                chrNum = loc[6]#[-1]
                if not pepdic[original_pep][1].has_key(chrNum):
                    pepdic[original_pep][1][chrNum] = [loc]
                else:
                    list_in_pepdic = pepdic[original_pep][1].get(chrNum)
                    is_in_list = 0
                    # for i in list_in_pepdic:
                    #     if i[3] == loc[3]:
                    #         is_in_list = 1
                    if is_in_list == 0:
                        pepdic[original_pep][1][chrNum].append(loc)

    del fasta_info
    del index_info
    del fasta_seq


    ####################################################################### grouping

    #location count fill in
    for pep in pepdic:
        count = 0
        for chr in pepdic[pep][1]:
            count += len(pepdic[pep][1][chr])
        if len(pepdic[pep]) == 5:
            pepdic[pep].append(count)
        else:
            pepdic[pep][5] = count


    #copy list
    coordi_list = {}
    for pep in pepdic:
        for chrNum in pepdic[pep][1]:
            if not chrNum in coordi_list.keys():
                coordi_list[chrNum] = {}
            for alist in pepdic[pep][1][chrNum]:
                coordi_list[chrNum][alist[3]] = pep

    pep_exons = defaultdict(lambda: defaultdict(list))

    for chrNum in coordi_list:
        list_in_chr = coordi_list[chrNum].keys()
        list_in_chr.sort()
        grouping_info = []
        gstart = -3000
        for index in range(len(list_in_chr)):
            if gstart + 5000 < list_in_chr[index]:
                grouping_info.append(index)
            gstart = list_in_chr[index]

    #     Sprob = []
    #     for group in range(len(grouping_info)-1):
    #         prob = 1
    #         for index in range(grouping_info[group],grouping_info[group+1]):
    #             prob *= (1- (1-pepdic[coordi_list[chrNum][list_in_chr[index]]][4]) / pepdic[coordi_list[chrNum][list_in_chr[index]]][5])
    #         prob = 1- prob
    #         Sprob.append(prob)
    #     prob = 1
    #     for index in range(grouping_info[len(grouping_info)-1],len(coordi_list[chrNum])):
    #         prob *= (1- (1-pepdic[coordi_list[chrNum][list_in_chr[index]]][4])/ pepdic[coordi_list[chrNum][list_in_chr[index]]][5])
    #     prob = 1- prob
    #     Sprob.append(prob)


        for coord in list_in_chr:
            Sp_index = binary_search(grouping_info,list_in_chr.index(coord),0,len(grouping_info)-1)
            pep = coordi_list[chrNum][coord]
            # output2.write(chrNum+'\t')
            for i in pepdic[pep][1].get(chrNum):
                if i[3] != coord:
                    continue
                for j in range(len(i[0])):
                    # output2.write(str(i[0][j])+'/'+str(i[1][j])+';')
                    pep_exons[pep][(chrNum,i[0][j],i[1][j])].append((i[7],i[8],i[9][j]))
                # output2.write('\t'+pep+'\t'+str(pepdic[pep][3])+'\t'+str(pepdic[pep][5])+'\t')
                # output2.write(str(pepdic[pep][4])+'\t'+str(i[-2])+'\n')
                #output.write(str(Sprob[Sp_index])+'\n')
    #             output.write(str(Sprob[Sp_index])+'\t'+str(i[-2])+'\n')

    return pep_exons

#functions
def binary_search(array,key,imin,imax):
    if (imax < imin):
        return imax
    else:
        imid = (imin+imax)/2
        if array[imid] > key:
            return binary_search(array,key,imin,imid-1)
        elif array[imid] < key:
            return binary_search(array,key,imid+1,imax)
        else:
            return imid


def    CleanPeptideString(pep_str):
    if pep_str[1] == '.':
        pep_str = pep_str[2:-2]
    pep_str = pep_str.replace("0","")
    pep_str = pep_str.replace("1","")
    pep_str = pep_str.replace("2","")
    pep_str = pep_str.replace("3","")
    pep_str = pep_str.replace("4","")
    pep_str = pep_str.replace("5","")
    pep_str = pep_str.replace("6","")
    pep_str = pep_str.replace("7","")
    pep_str = pep_str.replace("8","")
    pep_str = pep_str.replace("9","")
    pep_str = pep_str.replace("+","")
    pep_str = pep_str.replace(".","")
    pep_str = pep_str.replace("?","_")
    pep_str = pep_str.replace("_","")
    pep_str = pep_str.replace("(","")
    pep_str = pep_str.replace(")","")
    pep_str = pep_str.replace("[","")
    pep_str = pep_str.replace("]","")
    pep_str = pep_str.replace(",","")
    return pep_str

def get_location(start_seq,end_seq,start,end,length,strand,chrNum,gene,transcript):
    i = -1
    j = -1
    while start_seq >= 0:
        i += 1
        start_seq -= length[i]
    while end_seq > 0:
        j += 1
        end_seq -= length[j]
    if strand == 0:
        tmp_start = start[i] - start_seq
        tmp_end = start[j] - end_seq
    else:
        tmp_start = end[i] + start_seq
        tmp_end = end[j] + end_seq
    #print ' start,end,i,j ',tmp_start,tmp_end,i,j,
    location = [[],[],[]]
    if strand == 0:
        location.append(tmp_end)
        location.append(tmp_start)
        if j-i <= 0:
            location[0].append(tmp_end)
            location[1].append(tmp_start)
            location[2].append(tmp_start-tmp_end)
        else:
            location[0].append(start[i])
            location[1].append(tmp_start)
            location[2].append(tmp_start-start[i])
            for k in range(j-i-1):
                location[0].append(start[k+i+1])
                location[1].append(end[k+i+1])
                location[2].append(end[k+i+1]-start[k+i+1])
            location[0].append(tmp_end)
            location[1].append(end[j])
            location[2].append(end[j]-tmp_end)
    else:
        location.append(tmp_start)
        location.append(tmp_end)
        if j-i <= 0:
            location[0].append(tmp_start)
            location[1].append(tmp_end)
            location[2].append(tmp_end-tmp_start)
        else:
            location[0].append(tmp_start)
            location[1].append(end[i])
            location[2].append(end[i]-tmp_start)
            for k in range(j-i-1):
                location[0].append(start[k+i+1])
                location[1].append(end[k+i+1])
                location[2].append(end[k+i+1]-start[k+i+1])
            location[0].append(start[j])
            location[1].append(tmp_end)
            location[2].append(tmp_end-start[j])
    '''
    location = ''
    if strand == 0:
        if j-i <= 0:
            location = str(tmp_end) + '-' + str(tmp_start)
        else:
            location += str(start[i])+'-'+str(tmp_start)+';'
            for k in range(j-i-1):
                location += str(start[k+i]) +'-' + str(end[k+i]) +';'
            location += str(tmp_end) + '-'+str(end[j])
    else:
        if j-i <= 0:
            location = str(tmp_start) + '-' + str(tmp_end)
        else:
            location = str(tmp_start)+'-'+str(end[i])+';'
            for k in range(j-i-1):
                location += str(start[k+i]) +'-' + str(end[k+i]) +';'
            location += str(start[j]) + '-'+str(tmp_end)
    '''
    #pos 5
    location.append(strand)
    #pos 6
    location.append(chrNum)
    #pos 7
    location.append(gene)
    #pos 8
    location.append(transcript)
    #pos 9
    location.append([])
    for k in range(i,j+1):
        location[-1].append((start[k],end[k]))
    return location

sequences = []
with open(sys.argv[4]) as f:
    for l in f:
        sequences.append(l.rstrip())

exon_map = map_peptides_to_exons(sys.argv[1],sequences)

for peptide, mapped_exons in exon_map.items():
    print(peptide)
    for exon, mapping_locations in mapped_exons.items():
        print("\t{0}:{1}/{2}".format(*exon))
        for location in mapping_locations:
            print("\t\t{0} {1} ({2[0]}/{2[1]})".format(*location))
