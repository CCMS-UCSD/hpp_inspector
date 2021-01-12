'''
Created on Jul 7, 2017

@author: s3cha
'''
import re

class TagSearch(object):
    tag_table = {}
    tag_length = None
    def __init__(self):
        pass
    
    def makeTagTable(self,fasta,tag_length):
        self.tag_length = tag_length
        for index in range(len(fasta)-tag_length+1):
            tag = fasta[index:index+tag_length]
            if tag not in self.tag_table:
                self.tag_table[tag] = []
            self.tag_table[tag].append(index)
        pass
    
    def getAllLocation(self,peptide,fasta):
        if len(peptide) < self.tag_length:
            return [m.start() for m in re.finditer(peptide,fasta)]
        else:
            result = []
            tag = peptide[:self.tag_length]
            candidates = self.tag_table.get(tag)
            if candidates == None:
                return result
            for candi in candidates:
                if fasta[candi:candi+len(peptide)] == peptide:
                    result.append(candi)
            return result
                    
if __name__ == '__main__':
    x = TagSearch()
    known_DB = '/home/s3cha/data/Known_DB/Ensembl_GRCH38/Homo_sapiens.GRCh38.pep.all.fa'
    fasta_database = open(known_DB,'r')
    fasta_info = []
    index_info = [0]
    total_len = 0
    total_seq = []
    fseq = ''
    fhead = ''
    for line in fasta_database:
        if line.startswith('>'):
            if fseq != '':
                fasta_info.append(fhead)
                total_seq.append(fseq.replace('I','L'))
                total_len += len(fseq+'X')
                index_info.append(total_len)
            fhead = line.strip()
            fseq = ''
        else:
            fseq += line.strip()

    fasta_seq = 'X'.join(total_seq)
    
    tagS = TagSearch()
    tagS.makeTagTable(fasta_seq,4)
    
    result = tagS.getAllLocation('KDVNAAIATIKT'.replace('I','L'), fasta_seq) 
    print result, bool(result) 
        
    