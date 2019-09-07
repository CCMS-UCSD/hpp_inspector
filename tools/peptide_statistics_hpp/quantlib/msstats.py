from csv import DictWriter

def write_annotations(annotations, output_file):

    header = ['Raw.file','Condition','BioReplicate','Experiment','IsotopeLabelType']

    for a in annotations:
        print(a)
        a['Raw.file'] = a['Run']
        a['IsotopeLabelType'] = 'L'

    with open(output_file, 'w') as w:
        r = DictWriter(w, fieldnames = header, extrasaction='ignore')
        r.writeheader()
        for a in annotations:
            r.writerow(a)
