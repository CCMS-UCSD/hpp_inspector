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
    ids = {}
    if args.mztab != "" and args.mztab != None:
        for mztab_file in glob.glob(args.mztab + '/*'):
            ids = mztab.read(mztab_file, ids)
    peptide_to_psm = defaultdict(list)
    for filescan,psm in ids.items():
        peptide_to_psm[''.join([p for p in psm.sequence if p.isalpha()])].append((filescan,psm))
        
    with open(args.novel_coverage) as f, open(args.novel_psms) as fw:
        header = ['protein','filename','scan','sequence','charge','type']
        r = csv.DictReader(f, delimiter = '\t')
        w = csv.DictWriter(fw, delimiter = '\t', fieldnames = [''])
        for l in r:
            supporting_peptides = l['supporting_peptides'].split(',')
            novel_peptides = l['novel_peptides'].split(',')
            protein = l['protein']
            for peptide in supporting_peptides + novel_peptides:
                for (filescan,psm) in peptide_to_psm[peptide]:
                    w.write({
                        'protein':protein,
                        'filename':filescan[0],
                        'scan':filescan[1],
                        'sequence':psm.sequence,
                        'charge':psm.charge,
                        'type':'Supporting' if peptide in supporting_peptides else 'Novel'
                    })



if __name__ == '__main__':
    main()
