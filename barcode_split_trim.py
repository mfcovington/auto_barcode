#!/usr/bin/env python
import argparse
import fileinput
import re
import os
import sys
from pprint import pprint as pp

current_version = "v1.4.0"

def main():
    args = get_options()
    barcode_table = get_barcodes(args.list, args.barcode, args.id)
    # barcode_table = get_barcodes(0,'ACGTG','demo')
    barcode_length = validate_barcodes(barcode_table)
    directory, fq_name, barcode_name, unmatched_fq = open_fq_files(barcode_table, args.fastq, args.outdir, args.prefix, args.suffix, args.autoprefix, args.autosuffix, args.barcode, args.stats)
    test_write(barcode_table)
    split_trim_barcodes(args.fastq, barcode_table, barcode_length, args.notrim, args.stats)
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

def open_fq_files(barcode_table, fastq, outdir, prefix, suffix, autoprefix, autosuffix, barcode, stats):
    if outdir:
        directory = outdir
        filename = os.path.basename(fastq[0])
    else:
        directory, filename = os.path.split(fastq[0])

    if not os.path.exists(directory):
        os.makedirs(directory)

    if len(fastq) == 1:
        fq_name, fq_suffix = os.path.splitext(filename)
    else:
        fq_name = "multi_fq"

    if not stats:
        if autoprefix:
            prefix = fq_name + "."
        elif prefix:
            prefix += "."
        else:
            prefix = ""

        if suffix:
            suffix = "." + suffix
        else:
            suffix = ""

        barcode_name = os.path.basename(barcode);
        unmatched_fq = "{0}/unmatched.{1}fq_{2}.bar_{3}{4}.fq".format(directory, prefix, fq_name, barcode_name, suffix)
        open(unmatched_fq, 'w')

        for seq in dict.keys(barcode_table):
            if autosuffix:
                suffix = "." + seq
            barcode_id = barcode_table[seq]['id']
            fq_out = "{0}/{1}{2}{3}.fq".format(directory, prefix, barcode_id, suffix)
            barcode_table[seq]['fh'] = open(fq_out, 'w')

    return directory, fq_name, barcode_name, unmatched_fq;

def test_write(barcode_table):
    for seq in dict.keys(barcode_table):
        barcode_table[seq]['fh'].write('Here is some text for {0}'.format(seq))

def split_trim_barcodes(fastq, barcode_table, barcode_length, notrim, stats):
    total_matched = 0
    total_unmatched = 0
    barcodes_obs = {}
    good_read_id = re.compile('^@')
    fq_data = fileinput.input(fastq)
    for read_id in fq_data:
        seq = fq_data.readline()
        qual_id = fq_data.readline()
        qual = fq_data.readline()

        if not re.match(good_read_id, read_id):
            sys.exit("Encountered sequence ID ({0}) that doesn't start with '\@' on line {1} of FASTQ file: {2}...\nInvalid or corrupt FASTQ file?\n".format(read_id.rstrip('\n'), fq_data.filelineno() - 3, fq_data.filename()))

        if not len(seq) == len(qual):
            sys.exit("Encountered unequal sequence and quality lengths for read ({0}) starting on line {1} of FASTQ file: {2}...\nInvalid or corrupt FASTQ file?\n".format(read_id.rstrip('\n'), fq_data.filelineno() - 3, fq_data.filename()))


def close_fq_files(barcode_table):
    for seq in dict.keys(barcode_table):
        barcode_table[seq]['fh'].close()


if __name__ == '__main__':
    main()
