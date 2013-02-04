#!/usr/bin/env python

import sys
import argparse
import pprint
from collections import OrderedDict
import json
from sets import Set

pp = pprint.PrettyPrinter(indent=0)

parser = argparse.ArgumentParser(description='Description.')

parser.add_argument('gene_association_file')
parser.add_argument('idmapping_file')
parser.add_argument('go_slim_file')

header = ['DB', 'DB_Object_ID', 'DB_Object_Symbol', 'Qualifier', 'GO ID', 'DB:Reference', 'Evidence', 'With', 'Aspect', 'DB_Object_Name', 'Synonym', 'DB_Object_Type', 'Taxon_ID', 'Date', 'Assigned_By', 'Annotation_Extension', 'Gene_Product_Form_ID']

args = parser.parse_args()

mapping_dict = OrderedDict()

# go slim
go_slim_filter = OrderedDict()
f = open(args.go_slim_file, 'r')
allGoDict = {}

for line in f:
    term, extra = line.split('=>')
    term = term.strip()

    slims, ancestors = extra.split('//')

    if not term in go_slim_filter:
        go_slim_filter[term] = Set()

    for slim in slims.split():
        slim = slim.strip()
        go_slim_filter[term].add(slim)

# id mappings
f = open(args.idmapping_file, 'r')
for line in f:
    (id1, source, id2) = line.strip().split()
    mapping_dict[id1] = id2

# gene associaion
f = open(args.gene_association_file, 'r')
f.readline()

# print two lines for each gene, for enst and entrez
for line in f:
    line_sp = line.strip().split("\t")
    az = OrderedDict(zip(header, line_sp))
    az['-'] = '-'
    ls = [az[a] for a in ['DB_Object_ID', 'Taxon_ID', '-', 'GO ID', 'Evidence', 'Qualifier', '-', 'DB:Reference', 'Aspect']]
    ls.append(','.join(go_slim_filter.get(ls[3], ['None'])))

    # add enst mappings
    if ls[0] in mapping_dict:
        ls[0] = mapping_dict[ls[0]]
        print '\t'.join(ls)

    # add entrez mappings
    ls = [az[a] for a in ['DB_Object_Symbol', 'Taxon_ID', '-', 'GO ID', 'Evidence', 'Qualifier', '-', 'DB:Reference', 'Aspect']]
    ls.append(','.join(go_slim_filter.get(ls[3], ['None'])))
    print '\t'.join(ls)
