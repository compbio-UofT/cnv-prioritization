#!/usr/bin/env python

import sys
import argparse
import math
from collections import OrderedDict
import json
from random import choice
import random
import copy
import os
import cPickle as pickle

EXT_GZ = "gz"
SUFFIX_GZ = os.extsep + EXT_GZ
from contextlib import closing, contextmanager
from gzip import open as _gzip_open

def gzip_open(*args, **kwargs):
    return closing(_gzip_open(*args, **kwargs))
    # return _gzip_open(*args, **kwargs)

def maybe_gzip_open(filename, *args, **kwargs):
    if filename.endswith(SUFFIX_GZ):
        return gzip_open(filename, *args, **kwargs)
    else:
        return open(filename, *args, **kwargs)

parser = argparse.ArgumentParser(description='Annotate CNVs with similarity scores.')

parser.add_argument('weka_file', help="The arff file to annotate.")
parser.add_argument('sim_file', help="Similarity scores file.")
parser.add_argument('output', choices=['no_sim', 'by_sim'], help="Type of output: no_sim to not use similarity scores, by_sim otherwise.")
parser.add_argument('remove_no_feature_instances', choices=['remove', 'keep'], help="Set to remove no feature instances.")
parser.add_argument('neighbour_weight_function', help="The function to use to apply to the neighbour weights. e.g. [x, 2**x, x**2, sqrt(x)]")
parser.add_argument('weighted_gene_duplication', \
                    choices=['weighted_gene_duplication_no', 'weighted_gene_duplication_sim', 'weighted_gene_duplication_uniform', 'weighted_gene_duplication_max'], \
                    help="no to not use weights, sim to weight by similarity score, uniform for a uniform weighting out of 1000, max only use the highest similarity score")
parser.add_argument('plusone', choices=['plusone', 'pluszero'], \
                    help="plusone to use all instances even if the similarity score is 0, pluszero to not include these.")
parser.add_argument('out_gene_file', help="The output file, an arff annotated with similarity scores.")
parser.add_argument('similarity_rank_cutoff', type=int)
parser.add_argument('balance_test', choices=['bt_none', 'bt_remaining', 'bt_patient', 'bt_ptrem'], \
                    help="none for not balancing, remaining to balance and keep the remaining, patient to balance by patient, ptrem to balance by patient and keep the remainging.")
parser.add_argument('--balance_genes', '-b', action='store_true')
parser.add_argument('--debug', '-d', help="Debug.", action='store_true')

args = parser.parse_args()
if args.neighbour_weight_function == "0":
    anwf = "1"
else:
    anwf = args.neighbour_weight_function
fn_str = 'lambda x: %s' % args.neighbour_weight_function
print fn_str
NEIGHBOUR_WEIGHT_FUNCTION = eval(fn_str)

f = _gzip_open(args.weka_file, 'r')
line = f.readline().split('\t')
for i in range(len(line)):
    if '{' in line[i]:
        cutoff_index = i
        break
f.close()

prev = None
current = None

# out_gene_file = open(args.out_gene_file, 'w')

log_line_original_gene = []
log_score_original_gene = []
log_line_cnv = []
log_score_cnv = []

# out_cnv_file = open('weka_out_cnv.txt', 'w')
# discarded_cnvs_file = open('discarded_cnvs.txt', 'w')
# discarded_genes_file = open('discarded_genes.txt', 'w')

sim_file = _gzip_open(args.sim_file, 'r')

weka_file = _gzip_open(args.weka_file, 'r')
case_control = None

sim_line = sim_file.readline()
weka_line = weka_file.readline()

gene_to_sim = OrderedDict()

keep_looping = True

info_dict = OrderedDict()

while keep_looping:

    if not sim_line:
        keep_looping = False
        sim_line = '#'
        weka_line = '#'

    # ASSUME: line is not empty string
    # not original_gene line
    if sim_line[0] != "#":

        sim_line = sim_line.rstrip().split('\t')

        if args.debug:
            print 'sim_line:', sim_line

        gene_to_sim[sim_line[2]] = {'weight' : float(sim_line[3]),
                                    'rank' : int(sim_line[5]),
                                    'sim' : float(sim_line[12]),
                                    }

    # original_gene line
    else:

    	# ASSUME: line is not empty string
    	# not original_gene line
    	while weka_line[0] != '#':
            weka_line = weka_line.rstrip().split('\t')

            gene = weka_line[2]
            original_gene = weka_line[4]

            if args.debug:
                print '  each_gene_weka_line:', weka_line
                print '  gene:', gene
                print '  og:', original_gene

            # HARDCODE: phenotype similarity index
            if args.output == 'no_sim':
                sim_score = 0
            elif args.output == 'by_sim':
            	if args.debug:
                    print gene_to_sim
            	if gene_to_sim[gene]['rank'] > args.similarity_rank_cutoff:
                    sim_score = 0
                else:
                    sim_score = gene_to_sim[gene]['sim'] * NEIGHBOUR_WEIGHT_FUNCTION(gene_to_sim[gene]['weight'])
                # if gene in gene_to_sim: # TODO: need new similairty file!!
                #     sim_score = gene_to_sim[gene]['sim'] * NEIGHBOUR_WEIGHT_FUNCTION(gene_to_sim[gene]['weight'])
                # else:
                #     sim_score = 0
            replace_str = ', 6 %s' % sim_score

            out_line = '\t'.join(weka_line[cutoff_index:] + [str(sim_score)]).replace(',6 0', replace_str)

            # if specified in args, do not add no feature genes
            if not (args.remove_no_feature_instances == 'remove' and len(weka_line[0]) == 0):
                log_line_original_gene.append(out_line)
                log_score_original_gene.append(sim_score)
            
            weka_line = weka_file.readline()

        gene_to_sim = OrderedDict()

        # if not at the end of the file
        if weka_line != '#':
            # after while, original_gene_line
            weka_line = weka_line.rstrip().split('\t')
            # cnv name
            prev = current
            current = weka_line[9]

            if args.debug:
                print '  og_weka_line:', weka_line
                # print '    ', prev, current, current != prev
                # print '      log_line_original_gene', log_line_original_gene
                print '      log_score_original_gene', log_score_original_gene
        else:
            prev = current
            current = None

        # new cnv
        if current != prev and log_line_cnv != [] and case_control:

            # if case_control == 'BENIGN':
            #     log_score_cnv = [(max(log_score_cnv)-ls) for ls in log_score_cnv]

            max_score = max(log_score_cnv)
            idx = log_score_cnv.index(max_score)
            sum_sim_score = sum(log_score_cnv)

            if not case_control in info_dict[phenotype][patient_id]:
                info_dict[phenotype][patient_id][case_control] = OrderedDict()
            info_dict[phenotype][patient_id][case_control][prev] = []

            lsc_limit = float(1000)
            # lsc_limit = float(len(log_score_cnv) * 2)

            if ((sum_sim_score == 0 or len(log_score_cnv) > lsc_limit) and \
                (args.weighted_gene_duplication != 'weighted_gene_duplication_no')) or \
                ((args.weighted_gene_duplication == 'weighted_gene_duplication_no') and \
                 (case_control != 'HARMFUL' or args.output != 'by_sim')):
            	log_score_cnv_fractions = [1/lsc_limit for ls in log_score_cnv] # 1 per, special case (all 0 or > lsc_limit)
            elif args.weighted_gene_duplication == 'weighted_gene_duplication_uniform':
                log_score_cnv_fractions = [(1/float(len(log_score_cnv))) for ls in log_score_cnv] # uniform
            elif args.weighted_gene_duplication == 'weighted_gene_duplication_sim':
            	log_score_cnv_fractions = [ls/float(sum_sim_score) for ls in log_score_cnv] # fractions
            elif args.weighted_gene_duplication == 'weighted_gene_duplication_max' or \
                (args.weighted_gene_duplication == 'weighted_gene_duplication_no' and case_control == 'HARMFUL' and args.output == 'by_sim'):
            	log_score_cnv_fractions = [0 for ls in log_score_cnv]
                log_score_cnv_fractions[idx] = 1/lsc_limit
            else:
                raise Exception('Invalid args.weighted_gene_duplication. (%s, %s, %s)' % \
                                (args.weighted_gene_duplication, case_control, args.output))

            if args.debug:
                print '          case_control', case_control
                print '          idx', idx
                print '          sum_sim_score', sum_sim_score
                print '          log_score_cnv', log_score_cnv
                print '          log_score_cnv_fractions', log_score_cnv_fractions

            log_score_cnv_fractions = [int(math.ceil(ls*lsc_limit)) for ls in log_score_cnv_fractions]

            if args.debug:
                print '          log_score_cnv_fractions', log_score_cnv_fractions

            # add 1 to all so genes with hposim of 0 are still represented, +1, plusone
            if args.plusone == "plusone":
                log_score_cnv_fractions = [ls+1 for ls in log_score_cnv_fractions]
                if args.debug:
                    print '          log_score_cnv_fractions', log_score_cnv_fractions

            if args.debug:
                print '          log_score_cnv_fractions_sum', sum(log_score_cnv_fractions)

            # # gene creation
            # if args.weighted_gene_duplication == 'weighted_gene_duplication_no':
            #     # output all
            #     if case_control == 'BENIGN' or (case_control == 'HARMFUL' and args.output == 'no_sim'):
            #         # cnv
            #         # print >> out_cnv_file, 'xxx'.join(log_line_cnv)
            #         # gene
            #         for log_line_cnv_line in log_line_cnv:
            #             # print >> out_gene_file, log_line_cnv_line

            #             # info_dict
            #             gene = log_line_cnv_line.split('\t')[7]
            #             info_dict[phenotype][patient_id][case_control][prev].append(log_line_cnv_line)

            #             if args.debug:
            #                 print 'patient_id', patient_id
            #                 print 'prev', prev
            #                 print 'gene', gene
            #                 print 'case_control', case_control
            #                 print 'log_line_cnv_line', log_line_cnv_line

            #     # output only top by sim score
            #     elif case_control == 'HARMFUL' and args.output == 'by_sim':
            #         # discard cnv/genes if no gene has a similarity score > 0
            #         if max(log_score_cnv) == 0:
            #             # cnv
            #             print >> discarded_cnvs_file, 'xxx'.join(log_line_cnv)
            #             # gene
            #             for log_line_cnv_line in log_line_cnv:
            #                 print >> discarded_genes_file, log_line_cnv_line
            #         else:
            #             # cnv
            #             # print >> out_cnv_file, log_line_cnv[idx]
            #             # gene
            #             # print >> out_gene_file, log_line_cnv[idx]

            #             # info_dict
            #             gene = log_line_cnv_line.split('\t')[7]

            #             if args.debug:
            #                 print 'patient_id', patient_id
            #                 print 'prev', prev
            #                 print 'gene', gene
            #                 print 'case_control', case_control
            #                 print 'log_line_cnv_line', log_line_cnv[idx]
                            
            #             info_dict[phenotype][patient_id][case_control][prev].append(log_line_cnv[idx])
            # elif args.weighted_gene_duplication == 'weighted_gene_duplication_sim' or \
            #     args.weighted_gene_duplication == 'weighted_gene_duplication_uniform':
            #     args.weighted_gene_duplication == 'weighted_gene_duplication_max':
            if True:
                max_counter = 0
                for llc, lscf in zip(log_line_cnv, log_score_cnv_fractions):
                    # out_cnv_file.write('%sxxx' % llc)

                    gene = llc.split('\t')[7]

                    if case_control == 'BENIGN':
                        is_max = 'benign'
                    elif case_control == 'HARMFUL' and max_counter == idx:
                        is_max = 'max'
                    else:
                        is_max = 'notmax'

                    # for iteration in range(int(lscf)):
                    #     # gene
                    #     # print >> out_gene_file, llc
                    #     # cnv
                    #     # out_cnv_file.write('%sxxx' % llc)
                    #     # info_dict

                    #     info_dict[phenotype][patient_id][case_control][prev].append('%s\t%s' % (llc, is_max))

                    weighted_weka_line = '%s\t%s' % (llc, is_max)
                    weighted_weka_line = weighted_weka_line.replace('}\t%', '}, {%s}\t%%' % lscf)
                    info_dict[phenotype][patient_id][case_control][prev].append(weighted_weka_line)
                    if args.debug:
                        print 'wwl: ', weighted_weka_line
                    max_counter += 1
                # out_cnv_file.write('\n')
            else:
                raise 'Invalid option for args.weighted_gene_duplication.'

            # reset cnv
            log_line_cnv = []
            log_score_cnv = []

        # if not at the end of the file
        if weka_line != '#':

            case_control = weka_line[10]
            phenotype = weka_line[1]
            if not phenotype in info_dict:
                info_dict[phenotype] = OrderedDict()
            patient_id = weka_line[8]
            if not patient_id in info_dict[phenotype]:
                info_dict[phenotype][patient_id] = OrderedDict()
                info_dict[phenotype][patient_id]['BENIGN'] = OrderedDict()
                info_dict[phenotype][patient_id]['HARMFUL'] = OrderedDict()

            # add the max sim score of neighbour gene as sim score for the original_gene
            if log_score_original_gene == []:
                max_score_original_gene = 0
            else:
                max_score_original_gene = max(log_score_original_gene)

            # TODO: put this in later, check downstream affects
            out_line = '\t'.join(weka_line[cutoff_index:] + [str(max_score_original_gene), weka_line[4].split(',')[0], weka_line[1]]).replace(',6 0', ',6 %s' % max_score_original_gene)

            if not (args.remove_no_feature_instances == 'remove' and len(weka_line[0]) == 1):
                log_line_cnv.append(out_line)
                log_score_cnv.append(max_score_original_gene)

            # reset original_gene
            log_line_original_gene = []
            log_score_original_gene = []

            weka_line = weka_file.readline()

    sim_line = sim_file.readline()

weka_file.close()
sim_file.close()
# out_cnv_file.close()
# out_gene_file.close()

out_gene_file = open(args.out_gene_file, 'w')
pickle.dump(info_dict, out_gene_file)
out_gene_file.close()

# print "info_dict_res"
# print json.dumps(info_dict, indent=4)
