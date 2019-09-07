from csv import DictReader
from quantlib import psm

def read(mod_decode_file, ids):
    with open(mod_decode_file) as f:
        r = DictReader(f, delimiter = '\t')
        for l in r:
            filename = l['ipSpecFile'].split('.')[0]
            scan = int(l['ipScanID'])
            peptide = l['curatedPept']
            charge = int(l['CS'])
            type = "MOD_DECODE:{}".format(l['Known'])
            modification = l['ModAnnotation']
            prev_id = ids[(filename,scan)][0]
            if l['Known'] == 'UNIMOD':
                ids[(filename,scan)] = [psm.PSM(peptide,charge,type,modification,prev_id.rt,prev_id.protein)]
    return ids
