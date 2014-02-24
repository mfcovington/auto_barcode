#!/usr/bin/env python
import argparse
import re
import sys
from pprint import pprint as pp

current_version = "v1.4.0"

def main():
    args = get_options()
    barcode_table = get_barcodes(args.list, args.barcode, args.id)
    # barcode_table = get_barcodes(0,'ACGTG','demo')
    barcode_length = validate_barcodes(barcode_table)
    open_fq_files(barcode_table)
    test_write(barcode_table)
    # split_trim_barcodes()
    close_fq_files(barcode_table)
    # summarize_observed_barcodes()
    # summarize_counts()
    # plot_summary()
    pp(barcode_table)

def get_options():
    description = "Extracts fastq reads for specified barcode(s) from one or multiple FASTQ files."
    epilog = """
        An output file in fastq format is written for each barcode to the directory
        containing IN.FASTQ, unless an output directory is specified.
        The default name of the output file is ID.fq. The output names can be
        customized using the Naming Options.
    """
    parser = argparse.ArgumentParser(description=description, epilog=epilog)
    parser.add_argument('--version', action='version', version="%(prog)s {0}".format(current_version))
    parser.add_argument('-i', '--id', help='Sample or Experiment ID', required=True)
    parser.add_argument('-b', '--barcode', help='Barcode or file w/ list of barcodes to extract', required=True)
    parser.add_argument('-l', '--list', help='Indicate --barcode is a list of barcodes in a file', action="store_true")
    parser.add_argument('-n', '--notrim', help='Split without trimming barcodes', action="store_true")
    parser.add_argument('-s', '--stats', help='Output summary stats only (w/o creating fastq files)', action="store_true")
    parser.add_argument('-o', '--outdir', help='Output file is saved in the specified directory (or same directory as IN.FASTQ, if --outdir is not used)')
    parser.add_argument("fastq", help="FASTQ file(s). Use wildcards ('*') to match multiple input FASTQ files.", nargs="+")
    naming=parser.add_argument_group('optional naming arguments')
    naming.add_argument('--autoprefix', help='Append FASTQ file name onto output', action="store_true")
    naming.add_argument('--autosuffix', help='Append barcode onto output', action="store_true")
    naming.add_argument('--prefix', help='Add custom prefix to output')
    naming.add_argument('--suffix', help='Add custom suffix to output')
    return parser.parse_args()

def get_barcodes(list, barcode, id):
    barcode_table = {}
    if list == True:
        with open(barcode, 'r') as barcode_list:
            for line in barcode_list:
                seq, sample_id = line.split()
                barcode_table[seq] = {'id': sample_id, 'count': 0}
    else:
        barcode_table[barcode] = [id, 0]
    return barcode_table

def validate_barcodes(barcode_table):
    seqs = dict.keys(barcode_table)
    bad_seq = re.compile('[^ACGT]', re.IGNORECASE)
    if any(re.match(bad_seq, s) for s in seqs):
        sys.exit('Invalid barcode found!')

    min_length = len(min(seqs, key=len))
    max_length = len(max(seqs, key=len))
    if min_length != max_length:
        sys.exit("Unexpected variation in barcode length (min={0}, max={1})".format(min_length, max_length))
    return max_length

def open_fq_files(barcode_table):
    for seq in dict.keys(barcode_table):
        barcode_table[seq]['fh'] = open('temp/{0}'.format(seq), 'w')

def test_write(barcode_table):
    for seq in dict.keys(barcode_table):
        barcode_table[seq]['fh'].write('Here is some text for {0}'.format(seq))

def close_fq_files(barcode_table):
    for seq in dict.keys(barcode_table):
        barcode_table[seq]['fh'].close()


if __name__ == '__main__':
    main()
