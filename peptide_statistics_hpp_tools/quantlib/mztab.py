import csv
from quantlib import psm

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
        new_sequence = []
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
        'UNIMOD:35':'+15.994915'
    }
    if 'UNIMOD' in modification:
        return convert.get(modification)
    else:
        return modification.split(":")[1]

def read(mztab_file, ids):
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
        while(nextline[0:3] == 'PRH' or nextline[0:3] == 'PRT' or nextline == '\n'):
            nextline = f.readline()
        headers = nextline.rstrip().split('\t')
        mztab_dict = csv.DictReader(f, fieldnames = headers, delimiter = '\t')
        for row in list(mztab_dict):
            parent_mass = row.get('opt_global_precursor_neutral_mass',0)
            protein = row['accession']
            peptide = peptide_string(row['sequence'],row['modifications'])
            source_file, index = parse_spectrum_ref(row['spectra_ref'],filenames)
            rt = row.get('opt_global_RTMean')
            if rt:
                rt = float(rt)
            ids[(source_file, index)] = [psm.PSM(peptide, int(row['charge']), 'MZTAB', row['modifications'], rt, protein, parent_mass)]
    return ids
