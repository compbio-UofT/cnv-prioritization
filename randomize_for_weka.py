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
import re
from gzip import open as _gzip_open
import doctest

def enum(cls):
    return [t for t in dir(cls) if (not t.startswith('__')) and (not t.endswith('__'))]

class Bset:
    bt_cnv, bt_patient = range(2)

class Rset:
    no, new, last = range(3)

parser = argparse.ArgumentParser(description='Description.')

parser.add_argument('in_gene_file', help="Arff body file of genes.")
parser.add_argument('iteration', type=int, help="Number of iterations to run.")
parser.add_argument('out_file', help="Randomized output arff file of genes.")
parser.add_argument('header_file', help="Arff header file.")
parser.add_argument('balance_test', choices=enum(Bset), help="bt_cnv for balanced, bt_patient for unbalanced.")
parser.add_argument('num_sets', type=int, help="Number of sets to split the data into [x-fold].")
parser.add_argument('add_remaining', choices=enum(Rset), \
                    help="no for do not use the leftover data for another set, new for creating another set for leftover data, last for adding the remaining data to the last set")
parser.add_argument('weighted_gene_duplication', choices=['sim', 'uniform'], \
                    help="sim for using the similarity scores to weigh the instances, uniform for uniform weighting")
parser.add_argument('--balance_genes', '-b', action='store_true', help="Set to balance the genes between benign and harmful.")
parser.add_argument('--debug', '-d', action='store_true', help="Set to output debugging infomation.")
args = parser.parse_args()

def expand_info_dict(info_dict):
    clone_dict = OrderedDict()
    info_dict_cp = copy.deepcopy(info_dict)
    for phenotype_key, phenotype_val in info_dict_cp.iteritems():
        clone_dict[phenotype_key] = OrderedDict()
        for patient_key, patient_val in phenotype_val.iteritems():
            clone_dict[phenotype_key][patient_key] = OrderedDict()
            for case_control_key, case_control_val in patient_val.iteritems():
                clone_dict[phenotype_key][patient_key][case_control_key] = OrderedDict()
                for cnv_key, cnv_val in case_control_val.iteritems():
                    clone_dict[phenotype_key][patient_key][case_control_key][cnv_key] = []
                    for gene in cnv_val:
                        if gene.count('{') > 1:
                            weight = int(gene.split('{')[2].split('}')[0])
                            gene = re.sub('}, {.*}\t%', '}, {}\t%', gene)
                        else:
                            weight = 1
                            
                        for i in range(weight):
                            clone_dict[phenotype_key][patient_key][case_control_key][cnv_key].append(gene)
    return clone_dict

def query_info_dict(info_dict):
    counts = OrderedDict()
    counts['phenotype'] = 0
    counts['patient'] = 0
    counts['cnv'] = OrderedDict()
    counts['gene'] = OrderedDict()
    counts['gene']['HARMFUL'] = 0
    counts['gene']['BENIGN'] = 0
    counts['cnv']['HARMFUL'] = 0
    counts['cnv']['BENIGN'] = 0

    for phenotype_key, phenotype_val in info_dict.iteritems():
        counts['phenotype'] += 1
        for patient_key, patient_val in phenotype_val.iteritems():
            counts['patient'] += 1
            for case_control_key, case_control_val in patient_val.iteritems():
                for cnv_key, cnv_val in case_control_val.iteritems():
                    counts['cnv'][case_control_key] += 1
                    counts['gene'][case_control_key] += len(cnv_val)

    return counts


def test_one_phenotype_info_dict(info_dict):
    for patient_key, patient_val in info_dict['HP:0001508'].iteritems():
        print patient_key
        for case_control_key, case_control_val in patient_val.iteritems():
            print case_control_key
            for cnv_key, cnv_val in case_control_val.iteritems():
                print cnv_key
                print len(cnv_val)

# number of cnvs for a phenotype for benign/harmful
def num_per_phenotype(info_dict, phenotype_id, case_control_id):
    return sum([len(v[case_control_id].values()) if case_control_id in v else 0 for k,v in info_dict[phenotype_id].iteritems()])

# number of cnvs for a phenotype for a patient for benign/harmful 
def num_per_patient(info_dict, phenotype_id, patient_id, case_control_id):
    if case_control_id in info_dict[phenotype_id][patient_id]:
	return sum([len(v) for v in info_dict[phenotype_id][patient_id][case_control_id].values()])
    return 0

# does a particualr patient have both a harmful and a benign cnv?
def single_patient_per_phenotype_w_at_least_1_benign_and_1_harmful(info_dict, phenotype_id, patient_id):
    return (num_per_patient(info_dict, phenotype_id, patient_id, 'BENIGN') > 0) and \
        (num_per_patient(info_dict, phenotype_id, patient_id, 'HARMFUL') > 0)

# return the number of patients with at least 1 benign and 1 harmful cnv
# def num_patient_per_phenotype_w_at_least_1_benign_and_1_harmful(info_dict, phenotype_id, at_least):
#     return [single_patient_per_phenotype_w_at_least_1_benign_and_1_harmful(info_dict, phenotype_id, patient_id) \
#                  for patient_id in info_dict[phenotype_id].keys()].count(True) >= at_least

# return the patients with at least 1 benign and 1 harmful cnv
def ids_patient_per_phenotype_w_at_least_1_benign_and_1_harmful(info_dict, phenotype_id):
    return [patient_id \
            for patient_id in info_dict[phenotype_id].keys() \
            if single_patient_per_phenotype_w_at_least_1_benign_and_1_harmful(info_dict, phenotype_id, patient_id)]

# add a patient to cnv_ret
def add_patient(cnv_ret, info_dict, phenotype_id, patient_id, train_test):

    benign_cnvs = info_dict[phenotype_id][patient_id]['BENIGN'].keys()
    harmful_cnvs = info_dict[phenotype_id][patient_id]['HARMFUL'].keys()

    ben_num = len(benign_cnvs)
    harm_num = len(harmful_cnvs)

    cnv_ret['BENIGN'][train_test] = OrderedDict()
    cnv_ret['HARMFUL'][train_test] = OrderedDict()
    for cnv_id in benign_cnvs[:ben_num]:
        cnv_ret['BENIGN'][train_test][cnv_id] = info_dict[phenotype_id][patient_id]['BENIGN'][cnv_id]
    for cnv_id in harmful_cnvs[:harm_num]:
        cnv_ret['HARMFUL'][train_test][cnv_id] = info_dict[phenotype_id][patient_id]['HARMFUL'][cnv_id]

    del info_dict[phenotype_id][patient_id]

    return cnv_ret, info_dict

# add a cnv to cnv_ret
def add_cnv(cnv_ret, info_dict, phenotype_id, case_control_id, total_set_id):

    patient_id = choice([patient \
                          for patient in info_dict[phenotype_id].keys() \
                          if num_per_patient(info_dict, phenotype_id, patient, case_control_id) > 0])

    # if ben/harm exists for this patient
    if case_control_id in info_dict[phenotype_id][patient_id]:
        # if there are cnvs
        if info_dict[phenotype_id][patient_id][case_control_id].keys() != []:
            cnv_ret[case_control_id][total_set_id] = OrderedDict()
            cnv_id = choice(info_dict[phenotype_id][patient_id][case_control_id].keys())
            cnv_ret[case_control_id][total_set_id][cnv_id] = info_dict[phenotype_id][patient_id][case_control_id][cnv_id]
            del info_dict[phenotype_id][patient_id][case_control_id][cnv_id]
            total_set_id += 1
        # if there are no more cnvs, remove this ben/harm
        if info_dict[phenotype_id][patient_id][case_control_id].keys() == []:
            del info_dict[phenotype_id][patient_id][case_control_id]

    # if there are no more ben/harm, remove this patient
    if info_dict[phenotype_id][patient_id].keys() == []:
        del info_dict[phenotype_id][patient_id]

    return cnv_ret, info_dict, total_set_id

# add remaining patients/cnvs to the test set
def add_remaining_fn(cnv_ret, info_dict, phenotype_id, remaining_set):
    for patient_id, patient_val in info_dict[phenotype_id].iteritems():
        for case_control_id, case_control_val in patient_val.iteritems():
            if not cnv_ret[case_control_id][remaining_set]:
                cnv_ret[case_control_id][remaining_set] = OrderedDict()
            for cnv_id, cnv_val in case_control_val.iteritems():
                cnv_ret[case_control_id][remaining_set][cnv_id] = cnv_val
    return cnv_ret, info_dict

def fill_out_set(num_sets, fill):
    return dict([[sk, copy.deepcopy(fill)] for sk in range(num_sets)])

# generate train/test sets to output
def create_out_genes(cnv_ret, balance_genes, num_sets, balance_cnvs):
    set_keys = range(num_sets)
    out_genes = fill_out_set(num_sets, [])

    # for each fold
    for set_key in set_keys:
        ben_harm_list = {'BENIGN':[], 'HARMFUL':[]}
        loop_values = []

        if cnv_ret['BENIGN'][set_key]:
            benign_cnvs = cnv_ret['BENIGN'][set_key].keys()
            ben_num = len(benign_cnvs)
        if cnv_ret['HARMFUL'][set_key]:
            harmful_cnvs = cnv_ret['HARMFUL'][set_key].keys()
            harm_num = len(harmful_cnvs)

        # balance cnvs by benign/harmful for training
        if balance_cnvs:

            if cnv_ret['BENIGN'][set_key] and cnv_ret['HARMFUL'][set_key]:
            	cnv_num = min(ben_num, harm_num)
                ben_num = cnv_num
                harm_num = cnv_num
                random.shuffle(benign_cnvs)
                random.shuffle(harmful_cnvs)
            else:
            	cnv_ret['BENIGN'][set_key] = None
                cnv_ret['HARMFUL'][set_key] = None

        if cnv_ret['BENIGN'][set_key]:
            if args.debug and balance_cnvs and set_key == 2:
                print 'len ben: %s' % len(benign_cnvs)

            for cnv_id in benign_cnvs[:ben_num]:
                ben_harm_list['BENIGN'].append(cnv_ret['BENIGN'][set_key][cnv_id])

        if cnv_ret['HARMFUL'][set_key]:
            for cnv_id in harmful_cnvs[:harm_num]:
                ben_harm_list['HARMFUL'].append(cnv_ret['HARMFUL'][set_key][cnv_id])

        # balance benign/harmful genes
        if balance_genes and cnv_ret['BENIGN'][set_key] and cnv_ret['HARMFUL'][set_key]:
            new_ben = []
            new_harm = []
            for benign_list, harm_list in zip(ben_harm_list['BENIGN'], ben_harm_list['HARMFUL']):
                gene_limit = min(len(benign_list), len(harm_list))
                random.shuffle(benign_list)
                random.shuffle(harm_list)

                new_ben += [benign_list[:gene_limit]]
                new_harm += [harm_list[:gene_limit]]
            ben_harm_list['BENIGN'] = new_ben
            ben_harm_list['HARMFUL'] = new_harm

        # for train/test add all benign/harmful genes to one list
        if cnv_ret['BENIGN'][set_key]:
            loop_values += ben_harm_list['BENIGN']
        if cnv_ret['HARMFUL'][set_key]:
            loop_values += ben_harm_list['HARMFUL']
        for gene_list in loop_values:
            out_genes[set_key] += gene_list

    return out_genes

def is_cnv_ret_full(cnv_ret):
    ''' dict -> bool
    Return if cnv_ret has been filled completely by cnvs.  If true, distribute to files; otherwise, these cnvs should be put into the remaining set.
    >>> is_cnv_ret_full({'BENIGN':{'0':None, '1': None}, 'HARMFUL':{'0':None, '1': None}})
    False
    >>> is_cnv_ret_full({'BENIGN':{'0':True, '1': True}, 'HARMFUL':{'0':True, '1': True}})
    True
    '''
    ret = True
    for case_control_key, case_control_val in cnv_ret.iteritems():
        # print '    ', case_control_key,
        for cnv_set_key, cnv_set_val in case_control_val.iteritems():
            # print cnv_set_key, cnv_set_val != None,
            if cnv_set_val == None:
                ret = False
    # print
    return ret

def random_info_dict(info_dict, balance_test, num_sets):
    while info_dict.keys() != []:
        # reset for each phenotype
        if args.add_remaining == 'new':
            num_sets_w_new = num_sets + 1
        else:
            num_sets_w_new = num_sets

        cnv_ret = {'BENIGN':
                   fill_out_set(num_sets_w_new, None),
                   'HARMFUL':
                   fill_out_set(num_sets_w_new, None),
                   }

        if args.debug:
            print 'start'
        total_case_control_id = 'BENIGN'
        total_set_id = 0

        # always balance by phenotype
        phenotype_id = choice(info_dict.keys())
        if args.debug:
            print '  phenotype: %s' % phenotype_id

        # loop through random patients until able to add a harm/benign for train/test
        while total_set_id < num_sets or total_case_control_id != 'HARMFUL':

            # balance by patients
            if balance_test == 'bt_patient':
                # break loop if no more patients for this phenotype
                p_ids = ids_patient_per_phenotype_w_at_least_1_benign_and_1_harmful(info_dict, phenotype_id)
                if len(p_ids) <= 0:
                    break

                # choose random id from the list of pations with at least 1 ben/harm cnv
                patient_id = choice(p_ids)

                # add train patient
                # TODO: always unbalanced cnvs, add to output
                cnv_ret, info_dict = add_patient(cnv_ret, info_dict, phenotype_id, patient_id, total_set_id)

                total_case_control_id = 'HARMFUL'
                total_set_id += 1

            # balance by cnvs
            elif balance_test == 'bt_cnv':
                # break loop if no more cnvs for this phenotype/case_control
                if num_per_phenotype(info_dict, phenotype_id, total_case_control_id) < 1:
                    break

                cnv_ret, info_dict, total_set_id = add_cnv(cnv_ret, info_dict, phenotype_id, total_case_control_id, total_set_id)
                if total_set_id >= num_sets and total_case_control_id == 'BENIGN':
                    total_set_id = 0
                    total_case_control_id = 'HARMFUL'
                if args.debug:
                    print '      end  : ', total_case_control_id, total_set_id
            else:
                raise Exception('Invalid choice for balance_test: %s' % balance_test)

        if args.debug:
            print '        end inner', total_case_control_id, total_set_id

        if args.add_remaining == 'last':
            remaining_set = num_sets - 1
        elif args.add_remaining == 'new':
            remaining_set = num_sets
        elif args.add_remaining == 'no':
            remaining_set = num_sets
        else:
            raise Exception('Invaid choice for remaining: %s' % args.add_remaining)        

        # if not full ->
        # optional: add remaining cnvs to test
        if total_set_id < num_sets or total_case_control_id != 'HARMFUL':
            if args.add_remaining != 'no':
                cnv_ret, info_dict = add_remaining_fn(cnv_ret, info_dict, phenotype_id, remaining_set)

            # transfer to remaining, do not loop over remaing_set
            # ben/harm
            for case_control_id, case_control_val in cnv_ret.iteritems():
                if args.add_remaining != 'no':
                    if not cnv_ret[case_control_id][remaining_set]:
                        cnv_ret[case_control_id][remaining_set] = OrderedDict()
                # each set train/test/etc...
                for cnv_set_id in range(remaining_set):
                    cnv_set_val = cnv_ret[case_control_id][cnv_set_id]
                    # each cnv
                    if cnv_set_val:
                        for cnv_id, cnv_val in cnv_set_val.iteritems():

                            # remove the cnvs added to remaining
                            cnv_ret[case_control_id][cnv_set_id] = None

            del info_dict[phenotype_id]
            # print 'del: %s, %s' % (phenotype_id)

        # if there are cnvs to add to train/test
        if (total_set_id >= num_sets and total_case_control_id == 'HARMFUL') or args.add_remaining != 'no':
            if args.debug:
                print '          start return'
            balance_genes = args.balance_genes

            out_genes_cnv_bal = create_out_genes(cnv_ret, balance_genes, num_sets_w_new, True)
            out_genes_cnv_nbal = create_out_genes(cnv_ret, balance_genes, num_sets_w_new, False)

	    return phenotype_id, out_genes_cnv_bal, out_genes_cnv_nbal

    return None, {}, {}

print args.balance_test

def random_out(iteration, info_dict, num_sets):

    if args.add_remaining == 'new':
        num_sets_w_new = num_sets + 1
    else:
        num_sets_w_new = num_sets

    out_dict_encoding = fill_out_set(num_sets_w_new, {'cnvbal':{}, 'cnvnbal':{}})
    
    iteration_counts = OrderedDict()
    iteration_file = {}
    for set_id in range(num_sets_w_new):
        iteration_counts[set_id] = OrderedDict()
        iteration_file[set_id] = OrderedDict()
        for bal_nbal in ['cnvbal', 'cnvnbal']:
            iteration_counts[set_id][bal_nbal] = OrderedDict()
            iteration_counts[set_id][bal_nbal]['PHENOTYPE'] = OrderedDict()
            iteration_counts[set_id][bal_nbal]['PHENOTYPE']['HARMFUL'] = set()
            iteration_counts[set_id][bal_nbal]['PHENOTYPE']['BENIGN'] = set()
            iteration_counts[set_id][bal_nbal]['PATIENT'] = OrderedDict()
            iteration_counts[set_id][bal_nbal]['PATIENT']['HARMFUL'] = set()
            iteration_counts[set_id][bal_nbal]['PATIENT']['BENIGN'] = set()
            iteration_counts[set_id][bal_nbal]['CNV'] = OrderedDict()
            iteration_counts[set_id][bal_nbal]['CNV']['HARMFUL'] = set()
            iteration_counts[set_id][bal_nbal]['CNV']['BENIGN'] = set()
            iteration_counts[set_id][bal_nbal]['GENE'] = OrderedDict()
            iteration_counts[set_id][bal_nbal]['GENE']['HARMFUL'] = set()
            iteration_counts[set_id][bal_nbal]['GENE']['BENIGN'] = set()
            iteration_file[set_id][bal_nbal] = _gzip_open('%s_%s_%s_%s.arff.gz' % (set_id, args.out_file, iteration, bal_nbal), 'w')

    clone_dict = expand_info_dict(info_dict)
    # clone_dict = copy.deepcopy(info_dict)

    phenotype, out_genes_cnv_bal, out_genes_cnv_nbal = random_info_dict(clone_dict, args.balance_test, num_sets)
    i = 0

    while phenotype:
        # print query_info_dict(clone_dict)
        for bal_nbal, out_genes in zip(['cnvbal', 'cnvnbal'], [out_genes_cnv_bal, out_genes_cnv_nbal]):
            for label, gene_list in out_genes.iteritems():
                weka_gene_i = 0
                for weka_gene in gene_list:
                    (patient, cnv, case_control, dup, hposim, original_gene, hpo_term, infomax) = weka_gene.split('\t')[2:10]
                    iteration_counts[label][bal_nbal]['PHENOTYPE'][case_control].add(hpo_term)
                    iteration_counts[label][bal_nbal]['PATIENT'][case_control].add(patient)
                    iteration_counts[label][bal_nbal]['CNV'][case_control].add(cnv)
                    iteration_counts[label][bal_nbal]['GENE'][case_control].add('%s_%s_%s' % (cnv, original_gene, weka_gene_i))
                    # print hpo_term, patient, cnv, original_gene
                    if not weka_gene in out_dict_encoding[label][bal_nbal]:
                        out_dict_encoding[label][bal_nbal][weka_gene] = 0
                    out_dict_encoding[label][bal_nbal][weka_gene] += 1
                    weka_gene_i += 1


        phenotype, out_genes_cnv_bal, out_genes_cnv_nbal = random_info_dict(clone_dict, args.balance_test, num_sets)
        # print phenotype, query_info_dict(clone_dict)
        i += 1

    for bal_nbal in ['cnvbal', 'cnvnbal']:
        for label in range(num_sets_w_new):
            # num_count -> arff weight

            # TODO: maybe should be at randomization level?
            for gene_line, num_count in out_dict_encoding[label][bal_nbal].iteritems():
                if args.weighted_gene_duplication == 'sim':
                    exponator = 2
                else:
                    exponator = 1

                print >> iteration_file[label][bal_nbal], gene_line.replace('{}', '{%s}' % num_count**exponator)
                # TODO: use this to duplicate lines instead of weighing
                # for i in range(num_count):
                #     print >> iteration_file[label][bal_nbal], gene_line.replace('{}', '{%s}' % 1)
            iteration_file[label][bal_nbal].close()

    print 'iteration: %s, %s' % (iteration, i)
    for train_test_key, train_test_val in iteration_counts.iteritems():
        for subset_key, subset_val in train_test_val.iteritems():
            for case_control_key, case_control_val in subset_val.iteritems():
                for asdf_key, asdf_val in case_control_val.iteritems():
                    print train_test_key, subset_key, case_control_key, asdf_key, len(asdf_val)
    sys.stdout.flush()


in_gene_file = open(args.in_gene_file, 'r')
info_dict = pickle.loads(in_gene_file.read())
in_gene_file.close()

for iteration in range(1, args.iteration + 1):
    random_out(iteration, info_dict, args.num_sets)

doctest.testmod()
