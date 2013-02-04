#!/usr/bin/env python

import sys
import argparse
import os
import json
from collections import OrderedDict
import errno
from gzip import open as _gzip_open

def arg_parser():
    '''Parse the arguments of the script.
    '''
    parser = argparse.ArgumentParser(description='Find the CNVs which were originally labed as not harmful that are now classified as harmful.')
    parser.add_argument('all_db_file', help="Gff file of all CNVs.")
    parser.add_argument('out_pred_w_cnv_file', help="Weka predictions file with CNV annotations.")
    parser.add_argument('--debug', '-d', help="Debug.", action='store_true')

    args = parser.parse_args()
    return args

def load_cnv_dict(all_db_file):

    cnv_dict = OrderedDict()

    f = open(all_db_file, 'r')
    for line in f:
        line = line.rstrip()
        cnv_name = line.split()[2]
        cnv_dict[cnv_name] = line
        
    f.close()

    return cnv_dict

def find_harm_cnvs_in_no_harm_patients(out_pred_w_cnv_file, cnv_dict):
    f = _gzip_open(out_pred_w_cnv_file, 'r')
    for line in f:
        line = line.rstrip()
        line_split = line.split()
        nw_numharm = int(line_split[18])
        cnv_name = line_split[6]
        predicted = line_split[2]
        cnv_res = cnv_dict.get(cnv_name, [])
        cnv_res_split = cnv_res.split()
        if len(cnv_res) > 0:
            chrom = cnv_res_split[0]
            start = cnv_res_split[3]
            end = cnv_res_split[4]
        else:
	    chrom = 'None'
            start = 'None'
            end = 'None'

        if nw_numharm == 0 and predicted == '2:HARMFUL':
            print '\t'.join(['%s:%s-%s' % (chrom, start, end), line, cnv_res])
    f.close()
    
def main():
    '''Main.
    '''

    # get python params
    args = arg_parser()

    cnv_dict = load_cnv_dict(args.all_db_file)
    find_harm_cnvs_in_no_harm_patients(args.out_pred_w_cnv_file, cnv_dict)

main()
