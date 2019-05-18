import argparse
import sys
import csv
import glob
import mztab

def arguments():
    parser = argparse.ArgumentParser(description='mzTab to list of peptides')
    parser.add_argument('-m','--mztab', type = str, help='mzTab')
    parser.add_argument('-s','--spec_on_server', type = str, help='Spectra')
    parser.add_argument('-n','--novel_coverage', type = str, help='Novel Coverage')
    parser.add_argument('-p','--novel_psms', type = str, help='Novel PSMs')
    if len(sys.argv) < 4:
        parser.print_help()
        sys.exit(1)
    return parser.parse_args()

def main():
    args = arguments()


if __name__ == '__main__':
    main()
