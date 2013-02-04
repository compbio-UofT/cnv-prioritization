#!/usr/bin/env python

import sys
import argparse
import re
from gzip import open as _gzip_open

parser = argparse.ArgumentParser(description='Remove certain features from an arff file.')

parser.add_argument('--remove', '-R', help="The feature indexes to remove.")
parser.add_argument('--input', '-i', help="Input arff file.", required=True)
parser.add_argument('--output', '-o', help="Output arff file. If none, use stdout.")
parser.add_argument('--debug', '-d', help="Debug", action='store_true')

args = parser.parse_args()

if args.output == None:
    output_file = sys.stdout
else:
    output_file = _gzip_open(args.output, 'w')

intervals = args.remove.split(',')

remove_list = []
for interval in intervals:
    if '-' in interval:
        (start, end) = interval.split('-')
        start = int(start)
        end = int(end)
        remove_list += range(start-1, end)
    else:
    	remove_list.append(int(interval)-1)

remove_list_mone = [remove + 1 for remove in remove_list]

if args.debug:
    print remove_list

f = _gzip_open(args.input, 'r')
i = 0 # old index
j = 0 # new index
translation_mapping = {}
default_values = {}
for line in f:
    line = line.strip()

    # first line
    if '@relation' in line:
        print >> output_file, "%s -weka.filters.unsupervised.attribute.Remove-R%s" % (line, args.remove)
    # header
    elif '@attribute' in line:
        # arff header
        if not i in remove_list:
            print >> output_file, "%s\t%%\t%s\t%s" % (line, j, i)
            translation_mapping[i] = j
            feature_type = line.split()[2]
            if feature_type == 'NUMERIC':
                default_value = '0.0'
            elif '{' in feature_type and '}' in feature_type:
                default_value = feature_type[1:].split(',')[0]
            else:
                default_value = '0.0'

            default_values[i] = default_value
            j += 1
        i += 1

    # arff body
    elif len(line) > 0 and line[0] == '{':

        data_line = line.split('\t')[0]
        if data_line.count('{') > 1:
            split_idx = data_line.rfind(',')
            p1 = data_line[:split_idx]
            p2 = data_line[split_idx:]
            data_line = p1
            weighted = True
        else:
            weighted = False

        pair_list = []

        # shift arff feature indicies
        for pair in data_line[1:-1].split(','):
            (k, v) = pair.split()
            if int(k) in translation_mapping:
                # check values
                if args.debug:
                    print v, default_values[int(k)], v == default_values[int(k)]
                if not v == default_values[int(k)]:
                    pair_list.append(' '.join([str(translation_mapping[int(k)]), v]))
            # what is removed?
            else:
                if args.debug:
                    print k

        # arff body
        if weighted:
            print >> output_file, '{' + ','.join(pair_list) + '}%s\t%s' % (p2, '\t'.join(line.split('\t')[1:]))
        else:
            print >> output_file, '{' + ','.join(pair_list) + '}\t%s' % '\t'.join(line.split('\t')[1:])

    else:
        # non attribute/data lines
        print >> output_file, line
