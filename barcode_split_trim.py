#!/usr/bin/env python
from __future__ import division
import argparse
import fileinput
import re
import os
import numpy as np
import sys
from collections import Counter
from pandas import DataFrame
from rpy2 import robjects
from rpy2.robjects import pandas2ri
from rpy2.robjects.lib import ggplot2
from rpy2.robjects.packages import importr
from tabulate import tabulate
grdevices = importr('grDevices')
pandas2ri.activate()

current_version = "v1.4.0"

# Just a quick usage example:
# ./barcode_split_trim.py --id demo --list -b sample_files/barcode.list --out demo-output sample_files/sequences.fq

# TODO: Reformat per python standards
# TODO: Add NoPlot option (and put imports in the plot function?)

def main():
    args = get_options()
    barcode_table = get_barcodes(args.list, args.barcode, args.id)
    barcode_length = validate_barcodes(barcode_table)
    directory, fq_name, barcode_name, unmatched_fh = open_fq_files(barcode_table, args.fastq, args.outdir, args.prefix, args.suffix, args.autoprefix, args.autosuffix, args.barcode, args.stats)
    total_matched, total_unmatched, barcodes_obs = split_trim_barcodes(args.fastq, barcode_table, barcode_length, args.notrim, args.stats, unmatched_fh)
    close_fq_files(barcode_table, unmatched_fh)
    total_count = total_matched + total_unmatched
    summarize_observed_barcodes(barcode_table, barcodes_obs, total_count, directory, fq_name, barcode_name)
    summarize_counts(barcode_table, args.fastq, total_count, total_matched, total_unmatched, directory, fq_name, barcode_name)
    plot_summary(barcodes_obs, barcode_table, directory, args.id)

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

        barcode_name = os.path.basename(barcode)
        unmatched_fq = "{0}/unmatched.{1}fq_{2}.bar_{3}{4}.fq".format(directory, prefix, fq_name, barcode_name, suffix)
        unmatched_fh = open(unmatched_fq, 'w')

        for seq in dict.keys(barcode_table):
            if autosuffix:
                suffix = "." + seq
            barcode_id = barcode_table[seq]['id']
            fq_out = "{0}/{1}{2}{3}.fq".format(directory, prefix, barcode_id, suffix)
            barcode_table[seq]['fh'] = open(fq_out, 'w')

    return directory, fq_name, barcode_name, unmatched_fh

def split_trim_barcodes(fastq, barcode_table, barcode_length, notrim, stats, unmatched_fh):
    total_matched = 0
    total_unmatched = 0
    barcodes_obs = Counter()
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

        cur_barcode = seq[0:5]
        barcodes_obs[cur_barcode] += 1

        if barcode_table.has_key(cur_barcode):
            if not notrim:
                seq = seq[(barcode_length + 1):]
                qual = qual[(barcode_length + 1):]
            if not stats:
                barcode_table[cur_barcode]['fh'].write(read_id + seq + qual_id + qual)
            barcode_table[cur_barcode]['count'] += 1
            total_matched += 1
        else:
            if not stats:
                unmatched_fh.write(read_id + seq + qual_id + qual)
            total_unmatched += 1

    return total_matched, total_unmatched, barcodes_obs

def close_fq_files(barcode_table, unmatched_fh):
    unmatched_fh.close()
    for seq in dict.keys(barcode_table):
        barcode_table[seq]['fh'].close()

def summarize_observed_barcodes(barcode_table, barcodes_obs, total_count, directory, fq_name, barcode_name):
    rows = []
    for seq, count in barcodes_obs.most_common():
        sample = barcode_table.get(seq)['id'] if barcode_table.has_key(seq) else ""
        rows.append([seq, commify(count), percent(count, total_count), sample])
    table = tabulate(rows, headers=["barcode", "count", "percent", "id"], tablefmt="plain", stralign="right")
    summary = "{0}/log_barcodes_observed.fq_{1}.bar_{2}".format(directory, fq_name, barcode_name)
    with open(summary, 'w') as f:
        f.write(table)

def summarize_counts(barcode_table, fastq, total_count, total_matched, total_unmatched, directory, fq_name, barcode_name):
    summary = "{0}/log_barcode_counts.fq_{1}.bar_{2}".format(directory, fq_name, barcode_name)

    data = [('matched', total_matched), ('unmatched', total_unmatched)]
    table1 = table_id_num_pct(data, total_count)

    counts = map(lambda key: barcode_table[key]['count'], barcode_table.keys())
    data = [('min', min(counts)), ('max', max(counts)),
            ('mean', np.mean(counts)), ('median', np.median(counts))]
    table2 = table_id_num_pct(data, total_count)

    rows = []
    for seq in sorted(barcode_table, key=lambda key: barcode_table[key]['id']):
        sample = barcode_table[seq]['id']
        count = barcode_table[seq]['count']
        rows.append([sample, seq, commify(count), percent(count, total_count)])
    table3 = tabulate(rows, headers=["id", "barcode", "count", "percent"], tablefmt="plain", stralign="right")

    with open(summary, 'w') as f:
        f.write("Barcode splitting summary for:\n")
        for fq in fastq:
            f.write("  {0}\n".format(fq))
        f.write("---------------------------\n")
        f.write(table1)
        f.write("\n---------------------------\n")
        f.write("barcodes    {0}\n".format(len(barcode_table)))
        f.write(table2)
        f.write("\n---------------------------\n")
        f.write(table3)

def table_id_num_pct(data, total_count):
    rows = []
    for name, count in data:
        rows.append([name, commify(count), percent(count, total_count)])
    return tabulate(rows, tablefmt="plain", stralign="right")

def plot_summary(barcodes_obs, barcode_table, directory, expt_id):
    barcodes, counts, matches = get_vectors(barcodes_obs, barcode_table)
    df = DataFrame({'barcode': barcodes,
                    'count': counts,
                    'matched': matches})
    p = ggplot2.ggplot(df) + \
        ggplot2.aes_string(x='factor(matched)', y='count / 1000000') + \
        ggplot2.geom_boxplot(outlier_size = 0) + \
        ggplot2.geom_jitter() + \
        ggplot2.ggtitle(label = expt_id) + \
        ggplot2.ggplot2.xlab(label = "") + \
        ggplot2.scale_y_continuous(name = "Count\n(million reads)")

    filename = "{0}/{1}.png".format(directory, expt_id)
    grdevices.png(filename=filename, width=4, height=5, unit='in', res=300)
    p.plot()
    grdevices.dev_off()

def get_vectors(barcodes_obs, barcode_table):
    barcodes = []
    counts = []
    matches = []
    for barcode in dict.keys(barcodes_obs):
        barcodes.append(barcode)
        counts.append(barcodes_obs[barcode])
        match = "matched" if barcode_table.get(barcode) else "unmatched"
        matches.append(match)
    return barcodes, counts, matches

def percent(numerator, denominator, decimal_places=1):
    if denominator == 0:
        sys.exit("Oops, I divided by zero.")
    return "{0:.{1}f}%".format(100 * numerator / denominator, decimal_places)

def commify(number, decimal_places=0):
    return "{0:,.{1}f}".format(number, decimal_places)

if __name__ == '__main__':
    main()
