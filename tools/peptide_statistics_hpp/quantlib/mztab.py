import csv
from quantlib import psm
import math
import re

def filter_unknown_modification_rows(input_psms):
    output_psms = {}
    for scan, psm in input_psms.items():
        if not (psm[0].kind == 'MAESTRO' and psm[0].modifications != 'null'):
            output_psms[scan] = psm
    print("Removed {}/{} PSMs".format(len(input_psms)-len(output_psms), len(input_psms)))
    return output_psms

def peptide_string(sequence,modifications):
    if modifications == "null":
        return sequence
    else:
        new_sequence = sequence.split()
        mods = dict([
                (int(mod.split("-")[0]),find_mod(mod.split("-")[1]))
                for mod in modifications.split(",")
                if find_mod(mod.split("-")[1]) != None
            ])
        # mods = dict([
        #         (int(mod.split("-")[0]),find_mod(''.join(mod.split("-")[1:])))
        #         for mod in re.split("/\,(?![^\[]*\])/g", modifications)
        #         if find_mod(''.join(mod.split("-")[1:])) != None
        #     ])
        new_sequence = []
        new_sequence.append(mods.get(0,""))
        for (i,s) in enumerate(sequence):
            new_sequence.append(s)
            new_sequence.append(mods.get(i+1,""))
        return "".join(new_sequence)

def parse_spectrum_ref(spectrum_ref_string,filename_dict):
    fileref, index_string = spectrum_ref_string.split(":")
    index = int(index_string.replace("index=","").replace("scan=",""))
    filename = filename_dict.get(fileref)
    return (filename,index)

def find_mod(modification):
    convert = {
        'UNIMOD:1':'+42.010565',
        'UNIMOD:4':'+57.021464',
        'UNIMOD:5':'+43.005814',
        'UNIMOD:6':'+58.005479',
        'UNIMOD:7':'+0.984016',
        'UNIMOD:17':'+99.068414',
        'UNIMOD:21':'+79.966331',
        'UNIMOD:28':'-17.026549',
        'UNIMOD:34':'+14.015650',
        'UNIMOD:35':'+15.994915',
        'UNIMOD:214':'+144.102063',
        'UNIMOD:259':'+8.014199',
        'UNIMOD:267':'+10.008269',
        'UNIMOD:730':'+304.205360',
        'UNIMOD:731':'+304.199040',
        'UNIMOD:737':'+229.162932'
    }
    if 'UNIMOD' in modification:
        return convert.get(modification)
    elif 'PSI-MS' in modification:
        # psi_mod_split = modificatiom[1:-1].split(',')
        # if psi_mod_split[1] == 'MS:1001524':
        #     return '-{}'.format(psi_mod_split[3])
        # else:
            # print(modification)
        raise Exception
    else:
        # print(modification)
        return modification.split(":")[1]

def read(mztab_file, ids, mangled_name = None):
    filenames = {}
    with open(mztab_file) as f:
        nextline = f.readline()
        while(nextline[0:3] == 'MTD'):
            if 'ms_run' in nextline:
                header_line = nextline.rstrip().split("\t")
                ms_filename = header_line[1].replace("-location","")
                ms_filepath = header_line[2].replace("file://","f.")
                filenames[ms_filename] = ms_filepath
            nextline = f.readline()
        while(nextline[0:3] == 'PRH' or nextline[0:3] == 'PRT' or nextline[0:3] == 'COM' or nextline == '\n'):
            nextline = f.readline()
        headers = nextline.rstrip().split('\t')
        mztab_dict = csv.DictReader(f, fieldnames = headers, delimiter = '\t')
        for row in list(mztab_dict):
            parent_mass = row.get('opt_global_Precursor',0)
            protein = row.get('accession')
            try:
                peptide = peptide_string(row['sequence'],row['modifications'])
            except:
                print(row['sequence'],row['modifications'])
                raise Exception
            source_file, index = parse_spectrum_ref(row['spectra_ref'],filenames)
            rt = row.get('opt_global_RTMean')
            score = row.get('opt_global_EValue',0)
            if rt:
                rt = float(rt)
            if float(score) > 0:
                score = -math.log10(float(score))
            else:
                score = 100
            search_engine = 'MZTAB'
            # search_engine = row.get('search_engine','[,,MZTAB,]')[1:-1].split(',')[2]
            ids[(source_file, index)] = [psm.PSM(peptide, int(row['charge']), search_engine, row['modifications'], rt, protein, parent_mass, score, mangled_name)]
    return ids

def read_lib(mztab_file, ids):
    filenames = {}
    with open(mztab_file) as f:
        r = csv.DictReader(f, delimiter = '\t')
        for l in r:
            filename = l['filename']
            scan = int(l['scan'].replace('scan=',''))
            peptide = l['annotation']
            charge = int(l['charge'])
            parent_mass = float(l.get('mz',1))
            score = float(l['score'])
            ids[(filename, scan)] = [psm.PSM(peptide, charge, 'MSGF_AMB', ' ', None, None, parent_mass, score, filename)]
    return ids
