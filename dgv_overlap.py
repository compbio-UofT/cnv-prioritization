#!/usr/bin/env python

import sys
import argparse
from pandas import *
import itertools
import copy
from gzip import open as _gzip_open

parser = argparse.ArgumentParser(description='Calculate the DGV overlap feature for each CNV.')

parser.add_argument('cnvs_w_dgv_bed', help="Input CNVs/DGV overlap bed file.")
parser.add_argument('--debug', '-d', help="Debug.", action='store_true')

args = parser.parse_args()

# arff
f = _gzip_open(args.cnvs_w_dgv_bed, 'r')
prev_cnv = None
cnv_name = None
overlap_log = []

line = f.readline()
keep_looping = True

while keep_looping:

    if not line:
        keep_looping = False
        cnv_name = None
    else:
        line = line.strip().split('\t')
        prev_cnv = cnv_name
        cnv_name = line[2]
        extras = dict([pair.split('=') \
                        for pair in \
                        line[8].split(';')])

    # if at a new cnv, print the prev cnv
    if prev_cnv and prev_cnv != cnv_name:

        # take care of case where bedtools does not find an overlap (i.e. it puts dots at the end of the line)
        if sum(overlap_log) == 0:
            overlap_log = []

        # output
        print '\t'.join([str(s) \
                          for s in [prev_cnv, cnv_length, \
                                      len(overlap_log), \
                                      # METRIC: use the sum
                                      sum(overlap_log), \
                                      extras['SAMPLE_ID'], \
                                      extras['PHENOTYPE'][1:-1].split(',')[0].replace('CONTROL_', '')]]) # TODO: HARDCODE: assume 1 phenotype per cnv, take the first
        overlap_log = []

    # load up the next cnv
    if line:
        cnv_start = int(line[3])
        cnv_end = int(line[4])
        cnv_overlap = float(line[-1])
        cnv_length = cnv_end - cnv_start
        # METRIC: log the percentage overlap of the dgv cnv with our cnv
        if cnv_overlap/cnv_length >= .5:
            overlap_log.append(cnv_overlap/cnv_length)
        #     overlap_log.append(1)

        line = f.readline()

f.close()
