import argparse
import sys
from pathlib import Path
from csv import DictReader, DictWriter

def arguments():
    parser = argparse.ArgumentParser(description='mzTab to list of peptides')
    parser.add_argument('-i','--input_folder', type = Path, required=True, help='Input folder')
    parser.add_argument('-o','--output_file', type = Path, required=True, help='Output file')
    parser.add_argument('-d','--delimiter', type = str, help='Delimiter for file reading', default = '\t')

    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit(1)
    return parser.parse_args()

def main():

    args = arguments()

    all_fieldnames = set()

    for input_file in args.input_folder.glob('*'):
        with open(input_file) as input:
            r = DictReader(input, delimiter = args.delimiter)
            if r.fieldnames:
                all_fieldnames = all_fieldnames.union(set(r.fieldnames))

    with open(args.output_file, 'w') as output:
        w = DictWriter(output, delimiter = args.delimiter, fieldnames = list(all_fieldnames))
        w.writeheader()
        for input_file in args.input_folder.glob('*'):
            with open(input_file) as input:
                r = DictReader(input, delimiter = args.delimiter)
                for l in r:
                    w.writerow(l)


if __name__ == "__main__":
    main()
