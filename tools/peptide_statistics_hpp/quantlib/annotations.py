from csv import DictReader

def read(annotation_file, file_mapping):
    runs = []
    annotation_out = []
    with open(annotation_file) as f:
        r = DictReader(f)
        for l in r:
            annotation_out.append(l)
            mangled_name = None
            for key, value in file_mapping.items():
                if l['Run'].split('.')[0] in value:
                    mangled_name = key
            if mangled_name is None:
                print("File in annotations that isn't uploaded.")
            else:
                runs.append(l['Run'])
    return runs, annotation_out
