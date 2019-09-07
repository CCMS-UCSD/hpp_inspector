import xmltodict
from csv import DictReader
from collections import defaultdict

def read_params(input_file, full_paths = False):
    return get_mangled_file_mapping(parse_xml_file(input_file),full_paths)

def get_mangled_file_mapping(params, full_paths = False):
    all_mappings = params["upload_file_mapping"]
    mangled_mapping = {}
    for mapping in all_mappings:
        splits = mapping.split("|")
        if full_paths:
            mangled_name = splits[0]
            original_name = splits[1]
        else:
            mangled_name = splits[0].split(".")[0].split('/')[-1]
            original_name = splits[1].split(".")[0].split('/')[-1]
        mangled_mapping[mangled_name] = original_name

    return mangled_mapping

def parse_xml_file(input_file):
    with open(input_file) as f:
        key_value_pairs = defaultdict(list)
        xml_obj = xmltodict.parse(f.read())

        #print(json.dumps(xml_obj["parameters"]))
        for parameter in xml_obj["parameters"]["parameter"]:
            name = parameter["@name"]
            value = parameter["#text"]
            key_value_pairs[name].append(value)

    return key_value_pairs

def read_annotations(annotation_file, file_mapping):
    runs = []
    files_mangled = []
    with open(annotation_file) as f:
        r = DictReader(f)
        for l in r:
            mangled_name = None
            for key, value in file_mapping.items():
                if l['Run'].split('.')[0] in value:
                    mangled_name = key
            if mangled_name is None:
                print("File in annotations that isn't uploaded.")
            else:
                files_mangled.append(".".join(mangled_name.split("/")[-1].split(".")[:-1]))
                runs.append(l['Run'].split(".")[0].split('/')[-1])
    return runs, files_mangled
