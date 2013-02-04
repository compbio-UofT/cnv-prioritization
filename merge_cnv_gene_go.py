#!/usr/bin/env python

import sys
import csv
import cPickle as pickle
import networkx as nx
from sets import Set
import time
import operator
import argparse
from collections import OrderedDict
import copy
from operator import mul
import os

EXT_GZ = "gz"
SUFFIX_GZ = os.extsep + EXT_GZ
from contextlib import closing, contextmanager
from gzip import open as _gzip_open

def gzip_open(*args, **kwargs):
    return _gzip_open(*args, **kwargs)

def maybe_gzip_open(filename, *args, **kwargs):
    if filename.endswith(SUFFIX_GZ):
        return gzip_open(filename, *args, **kwargs)
    else:
        return open(filename, *args, **kwargs)

import gc
gc.disable()

parser = argparse.ArgumentParser(description='Formats the CNVs into an Arff. Anotates the CNVs/genes with go terms.')

parser.add_argument('gene_network_file', help="Gene-gene iteractions file.")
parser.add_argument('gene2go', help="Mapping between genes and go terms.")
parser.add_argument('case_control_gff', help="CNVs annotated with genes.")
parser.add_argument('weka_header_file', help="Output file. The weka header.")
parser.add_argument('go_slim_file', help="Mapping from go to go slim.")
parser.add_argument('evidence_code_filter_level', type=int, choices=[0,1,2,-1], \
                    help="0: no filter, n: tier filter, -1: modify weights")
parser.add_argument('filter_go_slim', choices=['y','n'], help="Toggle using go or go slim.")
parser.add_argument('take_top_x_neighbours', type=int, help=("Number of neighbours to use."))
parser.add_argument('neighbour_weight_function', help=\
                    "Function to modify the feature weights.")
parser.add_argument('go_and_ensg_to_hpo_file', help="Mapping between GO/genes/HPO.")
parser.add_argument('dijkstra_distance_cutoff', type=int, help="Max gene network distance.")
parser.add_argument('output_file', help="Output file. Weka body.")
parser.add_argument('--go_terms_for_sim', '-g', \
                    help="Output go terms for similarity processing.", action='store_true')

args = parser.parse_args()

totalstart = time.clock()
start = time.clock()

TIER_X_GO_FILTER = \
    OrderedDict([(-1, {}),
                  (0, {'EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'IEP', 'ISS', 'ISO', 'ISA', 'ISM'}),
                  (1, {'EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'IEP', 'ISS', 'ISO', 'ISA', 'ISM', \
                          'IGC', 'RCA', 'IBA', 'IBD', 'IKR', 'IRD', 'TAS', 'IC', 'IEA'})])

RESERVED_FIELD_NUM = 7
if args.neighbour_weight_function == "0":
    anwf = "1"
else:
    anwf = args.neighbour_weight_function
fn_str = 'lambda x: %s' % anwf
print fn_str
NEIGHBOUR_WEIGHT_FUNCTION = eval(fn_str)

output_file = _gzip_open(args.output_file, 'w')

def timer(s):
    global start
    global totalstart
    totalelapsed = (time.clock() - totalstart)
    elapsed = (time.clock() - start)
    start = time.clock()
    print '%s...%s...%s' % (s, elapsed, totalelapsed)
    sys.stdout.flush()

def load_gene_net():
    iin = open(args.gene_network_file, 'r')

    # new format
    # {'gene':[['neighbour':'weight'], ...]}
    gene_net = {}
    for line in iin:
        (gene, neighbours) = line.strip().split()
        neighbours = [(n.split(':')) for n in neighbours.split(',')]

        gene_net[gene] = neighbours[0:args.take_top_x_neighbours]

    iin.close()
    
    return gene_net

def load_gene2go():
    csvReader = csv.reader(open(args.gene2go, 'rb'), delimiter='\t')
    # build hash of refseq_ids -> go_terms
    gene2go = OrderedDict()
    for row in csvReader:
        gene = row[0]
        go = row[3]
        relation = row[5]

        if not 'NOT' in relation:
            if gene in gene2go:
                gene2go[gene].append(row)
            else:
                gene2go[gene] = [row]

    return gene2go

def load_go_slim_filter(RESERVED_FIELD_NUM):
    attribute_num = RESERVED_FIELD_NUM
    go_slim_filter = OrderedDict()
    f = open(args.go_slim_file, 'r')
    all_go_dict = {}

    for line in f:
        term, extra = line.split('=>')
        term = term.strip()
        
        slims, ancestors = extra.split('//')

        if not term in go_slim_filter:
            go_slim_filter[term] = Set()

        for slim in slims.split():
            slim = slim.strip()
            go_slim_filter[term].add(slim)

            # pre add to all_go_dict from slim file
            if not slim in all_go_dict and 'GO:' in slim:
                all_go_dict[slim] = attribute_num
                attribute_num += 1
        
    f.close()

    return go_slim_filter, all_go_dict, attribute_num

def load_ensg_to_hpo_file(f_name):
    ensg_to_hpo_dict = OrderedDict()
    f = open(f_name, 'r')
    for line in f:
        line = line.strip().split()
        enst = line[2]
        hpo_string = line[1]

        if enst != 'None':
           ensg_to_hpo_dict[enst] = Set(hpo_string.split(','))
    f.close()
    return ensg_to_hpo_dict

def add_go_and_parents(go_slim_filter, go_names, go_name, evidence_code_weight):
    if args.filter_go_slim == 'y':
        if go_name in go_slim_filter:
            for go_slim in go_slim_filter[go_name]:
                go_names[go_slim] = evidence_code_weight
    else:
        go_names[go_name] = evidence_code_weight

    return go_names

def inv_log(x):
    return 1/float(10**x)

def set_go_names(gene2go, gene_name, go_slim_filter):

    # build up list of go terms for this gene
    go_names = OrderedDict()
    for row in gene2go[gene_name]:
        go_name = row[3]
        evidenceCode = row[4]

        # if no evidence code filter level or if set to 0, do not filter
        if args.evidence_code_filter_level == 0:
            go_names = add_go_and_parents(go_slim_filter, go_names, go_name, 1)

        # modify weights
        elif args.evidence_code_filter_level == -1:
            for evidence_code_weight, evidence_code_tier in TIER_X_GO_FILTER.iteritems():
                if evidenceCode in evidence_code_tier:
                    go_names = add_go_and_parents(go_slim_filter, go_names, \
                                                 go_name, evidence_code_weight)
                    break

        # filter
        else:
            if evidenceCode in TIER_X_GO_FILTER.values()[args.evidence_code_filter_level]:
                go_names = add_go_and_parents(go_slim_filter, go_names, go_name, 1)

    return go_names

def update_genes_and_go(gene_name, go_names, all_go_dict, attribute_num, current_genes, \
                        original_gene, current_go_dict, neighbour_weight, neighbour_num):

    # for each go term
    for go_name, evidence_code_weight in go_names.iteritems():

        # increment all_go_dict if the current go term isn't already included
        if not go_name in all_go_dict:
            all_go_dict[go_name] = attribute_num
            attribute_num += 1

        code = all_go_dict[go_name]

        # increment current genes
        if not go_name in current_genes[original_gene][gene_name]['go']:
            current_genes[original_gene][gene_name]['go'][go_name] = {}
            current_genes[original_gene][gene_name]['go'][go_name][code] = []

        if args.neighbour_weight_function == "1":
            nwfnw = 1
        else:
            nwfnw = neighbour_weight
        current_genes[original_gene][gene_name]['go'][go_name][code].append\
            ([inv_log(neighbour_num), inv_log(evidence_code_weight), nwfnw])

        # increment current go dict
        # only add a go term once to a cnv-gene
        if not go_name in current_go_dict:
            current_go_dict[go_name] = 0
        current_go_dict[go_name] += evidence_code_weight

    return current_genes, current_go_dict, all_go_dict, attribute_num

def add_gene(go_slim_filter, gene2go, gene_name, current_genes, current_go_dict, all_go_dict, \
             attribute_num, neighbour_num, top_neighbour_num, neighbour_weight, \
             original_gene, current_gene_weights):

    current_gene_weights['weight'][gene_name] = neighbour_weight
    current_gene_weights['rank'][gene_name] = top_neighbour_num

    if not original_gene in current_genes:
        current_genes[original_gene] = {}

    # only add a gene once to a cnv
    if not gene_name in current_genes[original_gene]:
        current_genes[original_gene][gene_name] = {'go':{},
                                                 'weight':neighbour_weight,
                                                 'rank':top_neighbour_num,}

    # if the gene has some go terms
    if gene_name in gene2go:

        go_names = set_go_names(gene2go, gene_name, go_slim_filter)

        current_genes, current_go_dict, all_go_dict, attribute_num = \
            update_genes_and_go(gene_name, go_names, all_go_dict, attribute_num, current_genes, \
                        original_gene, current_go_dict, neighbour_weight, neighbour_num)

    return current_genes, current_go_dict, all_go_dict, attribute_num, current_gene_weights

def add_multiple(gene_name, gene_net, current_genes, gene2go, current_go_dict, all_go_dict, \
                 attribute_num, go_slim_filter, current_gene_weights):

    if gene_name in current_genes:
        # TODO:
        # print >> sys.stdout, "dup gene in this cnv root (%s)" % gene_name
        # do not allow multiple genes to be added with the same name within the same cnv
        # (gene dup, alt splicing, etc)
        return current_genes, current_go_dict, all_go_dict, attribute_num

    current_genes, current_go_dict, all_go_dict, attribute_num, current_gene_weights = \
        add_gene(go_slim_filter, gene2go, gene_name, current_genes, \
                 current_go_dict, all_go_dict, attribute_num, 0, 0, 1, gene_name, current_gene_weights)

    if args.take_top_x_neighbours >= 1:

        if gene_name in gene_net:

            top_neighbour_num = 1
            
            for neighbor_gene, neighbour_weight in gene_net[gene_name]:

                neighbour_weight = float(neighbour_weight)
                # HARDCODE: distance limit
                # TODO: breaks on dijkstra weight but sorted by RW weight!
                if neighbour_weight < 0.096/args.dijkstra_distance_cutoff:
                    break

                current_genes, current_go_dict, all_go_dict, attribute_num, current_gene_weights = \
                    add_gene(go_slim_filter, gene2go, neighbor_gene, \
                             current_genes, current_go_dict, all_go_dict, attribute_num, 1, \
                             top_neighbour_num, neighbour_weight, gene_name, current_gene_weights)

                top_neighbour_num += 1

    return current_genes, current_go_dict, all_go_dict, attribute_num

def set_printList(print_dict):
    sorted_print_dict = sorted(print_dict.iteritems(), key=operator.itemgetter(0))

    # create a list in sorted order
    printList = []
    for k, v in sorted_print_dict:

        # append to the sorted weka list
        printList.append("%s %s" % (k, v))

    return printList

def set_annotation_hpo_set(current_gene_list, ensg_to_hpo_dict):
    # TODO: "original" gene must not have a weight
    annotation_hpo_set = Set()
    for gene in current_gene_list:
        if gene in ensg_to_hpo_dict:
            annotation_hpo_set |= ensg_to_hpo_dict[gene]

    return annotation_hpo_set

def set_go_output_line(current_gene_list, gene2go):
    go_set = Set()
    for gene in current_gene_list:
        if gene in gene2go:
            for row in gene2go[gene]:
                go_name = row[3]

                # TODO: match evidence code level, do not use evidence code!
                # for modifying weights
                evidenceCode = row[4]

                if args.evidence_code_filter_level == 0:
                    go_set.add(go_name)
                else:
                    if evidenceCode in TIER_X_GO_FILTER.values()[args.evidence_code_filter_level]:
                        go_set.add(go_name)

    if args.go_terms_for_sim:
        go_output_line = ','.join(go_set)
    else:
        go_output_line = ','.join(go_set)[:2]

    return go_output_line

def output_line(print_dict, size, extras, col_dict, dup, phenotype_list, current_go_dict, \
                ensg_to_hpo_dict, info_lists, comment_prefix):

    current_gene_list = info_lists['gene_list']
    current_gene_list_weights = [str(w) for w in info_lists['gene_weight_list']]
    current_gene_list_ranks = [str(r) for r in info_lists['gene_rank_list']]
    original_genes = info_lists['original_gene_list']
    
    # sort the weka print list by attribute num
    if len(print_dict) > -1:

        printList = set_printList(print_dict)
        annotation_hpo_set = set_annotation_hpo_set(current_gene_list, ensg_to_hpo_dict)
        go_output_line = set_go_output_line(current_gene_list, gene2go)
            
        to_output_line = "%s%s\t%s\t%s\t%s\t%s\t%s\t{%s}\t%%\t%s\t%s\t%s\t%s" % \
            (comment_prefix, go_output_line, ','.join(phenotype_list), \
             ','.join(current_gene_list), ','.join(current_gene_list_weights), \
             ','.join(original_genes), ','.join(current_gene_list_ranks), \
             ",".join(printList), extras['SAMPLE_ID'], col_dict['NAME'], print_dict[0], dup)
        print >> output_file, to_output_line

def update_phenotypes(phenotypes):

    # for each phenotype of the current cnv
    phenotype_list = phenotypes[1:-1].split(",")

    # if there are more than one phenotype, mark as duplicated
    dup = False
    if len(phenotype_list) > 1:
        dup = True
    
    if 'CONTROL' in phenotypes or phenotypes == '[]':
        phenotype_list = [phenotype.strip()[8:] for phenotype in phenotype_list \
                          if phenotype != '']
        phenotype = 'BENIGN'
    else:
        phenotype_list = [phenotype.strip() for phenotype in phenotype_list if phenotype != '']
        phenotype = 'HARMFUL'

    return phenotype_list, dup, phenotype

def update_hpo(phenotype_list, all_hpo_dict, attribute_num, print_dict):
    for hpo in phenotype_list:

        hpo_list = [hpo]

        for hpo_query in hpo_list:

            if not hpo_query in all_hpo_dict:
                all_hpo_dict[hpo_query] = attribute_num
                attribute_num += 1

            # add total attributes it doesn't already exist
            # add to weka print list
            a_num = all_hpo_dict[hpo_query]
            print_dict[a_num] = 1

    return all_hpo_dict, attribute_num, print_dict

def update_info_lists(gene_id, original_gene, gene_list,
                      info_lists_cnv, info_lists_original_gene):
    info_lists_cnv['gene_list'].append(gene_id)
    info_lists_cnv['original_gene_list'].append(original_gene)
    info_lists_cnv['gene_weight_list'].append(gene_list['weight'])
    info_lists_cnv['gene_rank_list'].append(gene_list['rank'])

    info_lists_original_gene['gene_list'].append(gene_id)
    info_lists_original_gene['original_gene_list'].append(original_gene)
    info_lists_original_gene['gene_weight_list'].append(gene_list['weight'])
    info_lists_original_gene['gene_rank_list'].append(gene_list['rank'])

    info_lists_each_gene = {}
    info_lists_each_gene['gene_list'] = [gene_id]
    info_lists_each_gene['original_gene_list'] = [original_gene]
    info_lists_each_gene['gene_weight_list'] = [gene_list['weight']]
    info_lists_each_gene['gene_rank_list'] = [gene_list['rank']]

    return info_lists_cnv, info_lists_original_gene, info_lists_each_gene

def update_weights(go_id, go_count, all_go_dict, each_gene_dict, \
                   original_gene_dict, print_dict):
    # add total attributes it doesn't already exist
    # add to weka print list
    a_num = all_go_dict[go_id]

    go_tris = go_count[a_num]
    go_tris_mul = [reduce(mul, go_tri) for go_tri in go_tris]
    go_tris_mul_sum = sum(go_tris_mul)

    if args.neighbour_weight_function == "0":
        go_count_s = 1
    elif args.neighbour_weight_function == "1" or args.neighbour_weight_function == "x":
        go_count_s = go_tris_mul_sum
    else:
        go_count_s = NEIGHBOUR_WEIGHT_FUNCTION(go_tris_mul_sum)

    if not a_num in each_gene_dict:
        each_gene_dict[a_num] = 0
    each_gene_dict[a_num] += go_count_s

    if not a_num in original_gene_dict:
        original_gene_dict[a_num] = 0
    original_gene_dict[a_num] += go_count_s

    if not a_num in print_dict:
        print_dict[a_num] = 0
    print_dict[a_num] += go_count_s

    return each_gene_dict, original_gene_dict, print_dict, go_count_s

def start_current_genes_loop(base_dict, total_genes):

    info_lists_original_gene = {'gene_list':[],
                  'original_gene_list':[],
                  'gene_weight_list':[],
                  'gene_rank_list':[],}
    num_go_gene = 0
    sum_go_gene = 0

    original_gene_dict = copy.deepcopy(base_dict)
    total_genes += 1
    use_this_gene = False

    return info_lists_original_gene, num_go_gene, sum_go_gene, \
        original_gene_dict, total_genes, use_this_gene

def init_calc_weights(print_dict):

    base_dict = copy.deepcopy(print_dict)

    used_genes = 0
    total_genes = 0

    # add each go term to the weka list, get total position num from all_go_dict

    info_lists_cnv = {'gene_list':[],
                  'original_gene_list':[],
                  'gene_weight_list':[],
                  'gene_rank_list':[],}
    num_go_cnv = 0
    sum_go_cnv = 0

    return base_dict, used_genes, total_genes, info_lists_cnv, num_go_cnv, sum_go_cnv

def exit_current_genes_loop(use_this_gene, num_go_gene, sum_go_gene, \
                            used_genes, original_gene_dict):
    # TODO: use all genes, even those with no go features
    if use_this_gene:
        used_genes += 1

    # weight = current_gene_weights['weight'][gene_id]
    # HARDCODE
    original_gene_dict[4] = num_go_gene
    original_gene_dict[5] = sum_go_gene

    return used_genes, original_gene_dict

def current_genes_loop(original_gene_list, base_dict, use_this_gene, info_lists_cnv, \
                       info_lists_original_gene, original_gene, original_gene_dict, \
                       print_dict, num_go_cnv, num_go_gene, sum_go_cnv, sum_go_gene, \
                       size, extras, col_dict, dup, phenotype_list, \
                       current_go_dict, ensg_to_hpo_dict):

    for gene_id, gene_list in original_gene_list.iteritems():

        each_gene_dict = copy.deepcopy(base_dict)

        if gene_id == '.':
            continue
        else:
            use_this_gene = True

        info_lists_cnv, info_lists_original_gene, info_lists_each_gene = \
            update_info_lists(gene_id, original_gene, gene_list,
                              info_lists_cnv, info_lists_original_gene)

        for go_id, go_count in gene_list['go'].iteritems():

            each_gene_dict, original_gene_dict, print_dict, go_count_s = \
                update_weights(go_id, go_count, all_go_dict, each_gene_dict, \
                               original_gene_dict, print_dict)

            num_go_cnv += 1
            num_go_gene += 1
            sum_go_cnv += go_count_s
            sum_go_gene += go_count_s

        output_line(each_gene_dict, size, extras, col_dict, dup, phenotype_list, \
                    current_go_dict, ensg_to_hpo_dict, info_lists_each_gene, '')

    return use_this_gene, info_lists_original_gene, original_gene_dict, \
        num_go_cnv, num_go_gene, sum_go_cnv, sum_go_gene
    
def calc_weights(current_genes, print_dict, size, extras, col_dict, dup, phenotype_list, \
    current_go_dict, ensg_to_hpo_dict):

    base_dict, used_genes, total_genes, info_lists_cnv, num_go_cnv, sum_go_cnv = \
        init_calc_weights(print_dict)

    for original_gene, original_gene_list in current_genes.iteritems():

        info_lists_original_gene, num_go_gene, sum_go_gene, \
            original_gene_dict, total_genes, use_this_gene = \
            start_current_genes_loop(base_dict, total_genes)

        use_this_gene, info_lists_original_gene, original_gene_dict, \
            num_go_cnv, num_go_gene, sum_go_cnv, sum_go_gene = \
            current_genes_loop(original_gene_list, base_dict, use_this_gene, info_lists_cnv, \
                               info_lists_original_gene, original_gene, original_gene_dict, \
                               print_dict, num_go_cnv, num_go_gene, sum_go_cnv, sum_go_gene, \
                               size, extras, col_dict, dup, phenotype_list, \
                               current_go_dict, ensg_to_hpo_dict)

        used_genes, original_gene_dict = \
            exit_current_genes_loop(use_this_gene, num_go_gene, sum_go_gene, \
                                    used_genes, original_gene_dict)
        
        output_line(original_gene_dict, size, extras, col_dict, dup, phenotype_list, \
                    current_go_dict, ensg_to_hpo_dict, info_lists_original_gene, '#')
        
def process_cnv(col_dict, phenotype_set, current_go_dict, all_go_dict, current_genes, \
                all_hpo_dict, attribute_num, current_gene_weights, ensg_to_hpo_dict):

    extras = dict([extra.split('=') for extra in col_dict['GROUP'].split(';')])
    phenotypes = extras['PHENOTYPE']

    phenotype_list, dup, phenotype = update_phenotypes(phenotypes)

    # add the phenotype to the total set
    phenotype_set.add(phenotype)

    # start the list of attributes dict for weka
    size = int(col_dict['STOP']) - int(col_dict['START'])

    # TODO: this is counting duplicates now
    num_go = 0
    sum_go = 0

    # TODO: fix current_genes/current_go_dict!!
    # TODO: 3 doesn't make sense for genes!
    print_dict = {0 : phenotype, 1 : extras['TYPE'], 2 : size, \
                 3 : len(current_genes), 4 : num_go, 5 : sum_go, 6 : 0 }

    # add each hpo code to the weka list, get total position num from all_hpo_dict
    all_hpo_dict, attribute_num, print_dict = \
        update_hpo(phenotype_list, all_hpo_dict, attribute_num, print_dict)

    calc_weights(current_genes, print_dict, size, extras, col_dict, dup, phenotype_list, \
                 current_go_dict, ensg_to_hpo_dict)

    return attribute_num

def first_init():

    gff_header = ['CHR', 'SOURCE', 'NAME', 'START', 'STOP', 'SCORE', 'STRAND', \
                  'FRAME', 'GROUP', 'rCHR', 'rSTART', 'rSTOP', 'rGENE', 'rOVERLAP']
    
    phenotype_set = Set()
    last_col_dict = {}
    all_hpo_dict = {}

    return gff_header, phenotype_set, last_col_dict, all_hpo_dict

def initialize_vars():

    current_genes = OrderedDict()
    current_go_dict = {}
    current_gene_weights = OrderedDict()
    current_gene_weights['weight'] = OrderedDict()
    current_gene_weights['rank'] = OrderedDict()
    
    return current_genes, current_go_dict, current_gene_weights
        

def print_weka_body(gene_net, gene2go, go_slim_filter, all_go_dict, attribute_num, \
                    ensg_to_hpo_dict):
    
    # initialize
    gff_header, phenotype_set, last_col_dict, all_hpo_dict = first_init()
    current_genes, current_go_dict, current_gene_weights = initialize_vars()

    # parse csv
    csvReader = csv.reader(open(args.case_control_gff, 'rb'), delimiter='\t')
    for row in csvReader:
        col_dict = dict([(gff_header[i], row[i]) for i in range(len(row))])

        if 'NAME' in last_col_dict:

            if col_dict['NAME'] != last_col_dict['NAME']:
                attribute_num = process_cnv(last_col_dict, phenotype_set, \
                                            current_go_dict, all_go_dict, current_genes, \
                                            all_hpo_dict, attribute_num, current_gene_weights, \
                                            ensg_to_hpo_dict)

                current_genes, current_go_dict, current_gene_weights = initialize_vars()

        current_genes, current_go_dict, all_go_dict, attribute_num = \
            add_multiple(col_dict['rGENE'], gene_net, current_genes, gene2go, \
                         current_go_dict, all_go_dict, attribute_num, go_slim_filter, \
                         current_gene_weights)

        last_col_dict = col_dict

    attribute_num = process_cnv(last_col_dict, phenotype_set, current_go_dict, \
                                all_go_dict, current_genes, all_hpo_dict, attribute_num, \
                                current_gene_weights, ensg_to_hpo_dict)

    return all_go_dict, phenotype_set, all_hpo_dict

def print_weka_header(all_go_dict, phenotype_set, all_hpo_dict):
    header_file = open(args.weka_header_file, 'w')
    print >> header_file, "@relation 'cnvsk: -C -%d -split-percentage 50 \n" % len(all_hpo_dict)
    print >> header_file, "@attribute\tclass\t{%s}" % "BENIGN,HARMFUL"
    print >> header_file, "@attribute\tdupdel\t{%s}" % ",".join(['DUP','DEL'])
    print >> header_file, "@attribute\tsize\t%s" % "NUMERIC"
    print >> header_file, "@attribute\tnumgenes\t%s" % "NUMERIC"
    print >> header_file, "@attribute\tnumgo\t%s" % "NUMERIC"
    print >> header_file, "@attribute\tsumgo\t%s" % "NUMERIC"
    print >> header_file, "@attribute\thposimilarity\t%s" % "NUMERIC"
    merged_dict = dict(all_go_dict.items() + all_hpo_dict.items())
    sorted_merged_dict = sorted(merged_dict.iteritems(), key=operator.itemgetter(1))
    for k,v in sorted_merged_dict:
        if 'HP:' in k:
            print >> header_file, "@attribute\t%s\t{0,1}" % k
        else:
            print >> header_file, "@attribute\t%s\tNUMERIC" % k
    print >> header_file, ""
    print >> header_file, "@data"
    print >> header_file, ""
    header_file.close()


timer('Loading gene network')
gene_net = load_gene_net()

timer('Loading genes to GO')
gene2go = load_gene2go()

timer('Loading go slim filter')
go_slim_filter, all_go_dict, attribute_num = load_go_slim_filter(RESERVED_FIELD_NUM)

timer('Load ensg_to_hpo_dict')
ensg_to_hpo_dict = load_ensg_to_hpo_file(args.go_and_ensg_to_hpo_file)

timer('Print weka body')
all_go_dict, phenotype_set, all_hpo_dict = \
    print_weka_body(gene_net, gene2go, go_slim_filter, all_go_dict, \
                    attribute_num, ensg_to_hpo_dict)

timer('Print weka header')
print_weka_header(all_go_dict, phenotype_set, all_hpo_dict)

timer('Done')
