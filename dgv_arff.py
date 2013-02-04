#!/usr/bin/env python

import sys
import argparse
from pandas import *
import itertools
import copy
from gzip import open as _gzip_open
from collections import OrderedDict

parser = argparse.ArgumentParser(description="Create a CNV arff file from a gene arff file.")

parser.add_argument('cnvs_w_dgv_overlap', help="DGV annotation file.")
parser.add_argument('out_pt', help="Results pivot table file.  Summary of harmful/benign classification of genes per cnv.")
parser.add_argument('out_pred_w_cnv', help="Weka predictions with cnv annotations.")
parser.add_argument('--incorrect', '-i', help="Only output incorrect.", action='store_true')
parser.add_argument('--debug', '-d', help="Debug", action='store_true')

args = parser.parse_args()

cnv_dict = OrderedDict()
# arff
f = _gzip_open(args.cnvs_w_dgv_overlap, 'r')
for line in f:
    (cnv, length, num_overlap, metric_overlap, sample, phenotype) = line.rstrip().split()
    cnv_dict[cnv] = OrderedDict([['length', length],
                                  ['num_overlap', num_overlap],
                                  ['metric_overlap', metric_overlap],
                                  ['sample', sample],
                                  ['phenotype', phenotype],
                                  ['maxgene', 'NoCnv'],
                                  ['conf', -1],
                                  ['simscore', -1],
                                  ])
f.close()

# assign the max gene to look at for each cnv
f = _gzip_open(args.out_pred_w_cnv)
for line in f:
    line = line.rstrip().split()

    simscore = float(line[9])
    cnv = line[6]
    tmp = cnv_dict.get(cnv, OrderedDict())
    conf = float(line[4])

    if (simscore > tmp.get('simscore', -1)) or \
        ((simscore == tmp.get('simscore', -1)) and \
         (conf > tmp.get('conf', -1))):

        maxgene = line[10]
        tmp['maxgene'] = maxgene
        tmp['conf'] = conf
        tmp['simscore'] = simscore

        if cnv in cnv_dict:
            if args.incorrect and line[3] == '-':
                del cnv_dict[cnv]
            else:
                cnv_dict[cnv] = tmp
f.close()

print '@relation \'cnvsk: -C -23 -split-percentage 50'
print
print '@attribute\tclass\t{BENIGN,HARMFUL}'
print '@attribute\tlength\tNUMERIC'
print '@attribute\toverlap\tNUMERIC'
print '@attribute\tharmgene\tNUMERIC'
print
print '@data'
print

f = _gzip_open(args.out_pt, 'r')
f.readline()
f.readline()
f.readline()
for line in f:
    # actual_predicted
    (cnv, benign_benign, benign_harmful, harmful_benign, harmful_harmful) = line.rstrip().split()

    if 'Copy' in cnv:
        continue

    if int(benign_harmful) > 0 or int(harmful_harmful) > 0:
        harmgene = '1'
    else:
        harmgene = '0'

    if 'CNV' in cnv:
        case_control = 'HARMFUL'
        max_val = 'max'
    else:
        case_control = 'BENIGN'
        max_val = 'benign'

    if (not args.incorrect) or (cnv in cnv_dict):
        print '{0 %s, 1 %s, 2 %s, 3 %s}, {1}\t%%\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % \
            (case_control, cnv_dict[cnv]['length'], cnv_dict[cnv]['metric_overlap'], harmgene, \
             cnv_dict[cnv]['sample'], cnv, case_control, False, 0.0, None, cnv_dict[cnv]['phenotype'], max_val,
             cnv_dict[cnv]['maxgene'], cnv_dict[cnv]['conf'], cnv_dict[cnv]['simscore'])
f.close()
