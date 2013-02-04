#!/usr/bin/env python

import sys
import argparse
import pprint
from collections import OrderedDict
import json
from sets import Set

pp = pprint.PrettyPrinter(indent=0)

parser = argparse.ArgumentParser(description='Create an ENST to HPO mapping file.')

parser.add_argument('gene_to_hpo_file', help="genes_to_phenotype.txt from HPO.")
parser.add_argument('idmapping_file_ensg', help="idmapping of gene ids, filtered for ensg.")
parser.add_argument('idmapping_file_entrez', help="idmapping of gene ids, filtered for entrez.")
parser.add_argument('ensg_or_entrez_to_go_file', help="Mapping file between genes and go terms.")

args = parser.parse_args()

print >> sys.stderr, 'Loading gene to go'
if args.ensg_or_entrez_to_go_file != 'none':
    gene_to_go_set = {}
    f = open(args.ensg_or_entrez_to_go_file, 'r')
    # TODO: consider go evidence code?
    for line in f:
        line = line.strip().split('\t')
        gene = line[0]
        go = line[3]
        relation = line[5]

        if (not 'NOT' in relation) and (not 'ENSG' in gene):
            if not gene in gene_to_go_set:
                gene_to_go_set[gene] = Set()
            gene_to_go_set[gene].add(go)

print >> sys.stderr, 'Loading idmapping ensg'
mapping_dict_ensg = OrderedDict()
f = open(args.idmapping_file_ensg, 'r')
for line in f:
    (id1, source, id2) = line.strip().split('\t')
    mapping_dict_ensg[id1] = id2

print >> sys.stderr, 'Loading idmapping entrez'
mapping_dict_entrez = OrderedDict()
f = open(args.idmapping_file_entrez, 'r')
for line in f:
    (id1, source, id2) = line.strip().split('\t')
    mapping_dict_entrez[id2] = id1

f = open(args.gene_to_hpo_file, 'r')
f.readline()
go_to_hpo_dict = {}

print >> sys.stderr, 'Loop through annotations'
for line in f:
    (gene_id_entrez, gene_symbol, hpo_terms) = line.strip().split("\t")

    hpo_list = [hpo[-11:-1] for hpo in hpo_terms.split(';')]

    if args.ensg_or_entrez_to_go_file == 'none':

        # add enst mappings
        if gene_id_entrez in mapping_dict_entrez:
            gene_id_uniprot = mapping_dict_entrez[gene_id_entrez]
            if gene_id_uniprot in mapping_dict_ensg:
                gene_id_ensg = mapping_dict_ensg[gene_id_uniprot]

                print '\t'.join([gene_id_ensg, ','.join(hpo_list), gene_id_entrez, gene_symbol])
    else:
        if gene_symbol in gene_to_go_set:
            if gene_id_entrez in mapping_dict_entrez:
                gene_id_uniprot = mapping_dict_entrez[gene_id_entrez]
                if gene_id_uniprot in mapping_dict_ensg:
                    gene_id_ensg = mapping_dict_ensg[gene_id_uniprot]
                else:
                    gene_id_ensg = 'None'
            else:
                gene_id_ensg = 'None'
            print "%s\t%s\t%s\t%s" % (','.join(gene_to_go_set[gene_symbol]), ','.join(hpo_list), gene_id_ensg, gene_id_uniprot)
