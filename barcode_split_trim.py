#!/usr/bin/env python
from pprint import pprint as pp

def main():
    # validate_options()
    barcode_table = get_barcodes(1,'sample_files/barcode.list','demo')
    # barcode_table = get_barcodes(0,'ACGTG','demo')
    # barcode_length = validate_barcodes(barcode_table)
    # open_fq_fhs_in_bulk()
    # split_trim_barcodes()
    # close_fq_fhs()
    # summarize_observed_barcodes()
    # summarize_counts()
    # plot_summary()
    pp(barcode_table)

def get_barcodes(list, barcode, id):
    barcode_table = {}
    if list == True:
        barcode_list = open(barcode, 'r')
        for line in barcode_list:
            seq, sample_id = line.split()
            barcode_table[seq] = {'id': sample_id, 'count': 0}
        barcode_list.close()
    else:
        barcode_table[barcode] = [id, 0]
    return barcode_table

if __name__ == '__main__':
    main()
