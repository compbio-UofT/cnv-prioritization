#!/usr/bin/env python

import sys
import argparse
import networkx
import json
import os
import binascii
from gzip import open as _gzip_open

parser = argparse.ArgumentParser(description='Parse a j48 graph file into a cytoscape file.')

parser.add_argument('j48_graph_file', help="Input j48 graph file. Use '_' to create the namespace_dict file.")
parser.add_argument('go_file', help="GO file.")
parser.add_argument('hp_file', help="HPO file.")
parser.add_argument('--debug', '-d', help="Debug.", action='store_true')

args = parser.parse_args()

G = networkx.DiGraph()

uniq_id = binascii.hexlify(os.urandom(16))

# method to parse ontology files to dictionaries
# key: ontology terms
# value: list of ontology term parents
def parse_ontology(ont_file):
    ont_dict = {}
    namespace_dict = {}
    f = open(ont_file, 'r')
    current_ont = None
    current_name = None
    current_namespace = None
    for line in f:
        line = line.strip().split()
        if current_ont == 'Next':
            current_ont = line[1]
            # if not current_ont in ont_dict:
            #     ont_dict[current_ont] = []
        else:
            if len(line) >= 1:
                if line[0] == '[Term]':
                    current_ont = 'Next'
                elif len(line) >= 2:
                    if 'name:' == line[0]:
                        current_name = " ".join(line[1:])
                        ont_dict[current_ont] = current_name
                        # print "%s -> %s" % (current_ont, val)
                    elif 'alt_id:' == line[0]:
                        ont_dict[line[1]] = current_name
                    	namespace_dict[current_ont] = current_namespace
                        # print "%s -> %s" % (current_ont, val)
                    elif 'namespace:' == line[0]:
                        current_namespace = " ".join(line[1:])
                    	namespace_dict[current_ont] = current_namespace
    f.close()
    return ont_dict, namespace_dict

# file to create hp <-> hp layer
hp2parents, _ = parse_ontology(args.hp_file)

# file to create go <-> go layer
go2parents, namespace_dict = parse_ontology(args.go_file)

if args.j48_graph_file == '_':
    print json.dumps(namespace_dict, indent=4)

if not args.j48_graph_file == '_':

    f = _gzip_open(args.j48_graph_file, 'r')
    for line in f:
        line = line.strip()

        if not len(line) > 0:
            continue
        if not line[0] == 'N':
            continue

        if '->' in line:
            n1 = uniq_id + "_" + line.split('->')[0]
            n2 = uniq_id + "_" + line.split('->')[1].split()[0]
            label = line.split('"')[1]

            label = label.replace('>', 'gt')
            label = label.replace('<', 'lt')

            G.add_edge(n1, n2, label=label)
            # print 'edge: (%s, %s) [%s]' % (n1, n2, label)
        else:
            node_name = uniq_id + "_" + line.split()[0]
            node_label = line.split('"')[1]

            print_label = node_label

            if node_label in go2parents:
                print_label += " %s" % go2parents[node_label]
            if node_label in hp2parents:
                print_label += " %s" % hp2parents[node_label]
            if node_label in namespace_dict:
                print_label += " %s" % namespace_dict[node_label]

            G.add_node(node_name, label=print_label)
            # print 'node: (%s) [%s]' % (node_name, node_label)

    networkx.write_graphml(G, 'j48.xml')
