#!/usr/bin/env python

import sys
import argparse
from pandas import *
import itertools
import copy
from gzip import open as _gzip_open
from collections import OrderedDict

parser = argparse.ArgumentParser(description='Calculate correct, precision, recall, etc on results.')

parser.add_argument('test_arff_file', help="Arff file.")
parser.add_argument('pred_file', help="Weka predictions on the arff file.")
parser.add_argument('out_file_postfix', help="Postfix of output files.")
parser.add_argument('--debug', '-d', help="Debug.", action='store_true')

args = parser.parse_args()

patient_harm_ben_counts = OrderedDict()

# arff
arff_info_pre = [None]
f = _gzip_open(args.test_arff_file, 'r')
for preline in f:
    line = preline.strip().split('\t')

    # for confusion matrix
    if '@attribute' in line and 'class' in line:
        class_values = line[2][1:-1].split(',')
        class_values = ['%s:%s' % (i, val) for i, val in zip(range(1, len(class_values) + 1), class_values)]

    # create matrix of info in arff file
    elif len(line) > 0 and len(line[0]) > 0 and line[0][0] == '{':

    	if True:
    	# if not 'cnv' in args.test_arff_file:
            if preline.count('{') > 1:
                weight = int(preline.split('{')[2].split('}')[0])
            else:
            	weight = 1

            # HARDCODE!
            if len(line) == 10:
                arff_info_pre.append(line[2:10] + ['None', -1, -1] + [weight])
            elif len(line) == 13:
                arff_info_pre.append(line[2:13] + [weight])
            else:
            	raise(Exception("Wrong number of lines in arff file: %s" % len(line)))

            patient = line[2]
            case_control = line[4]

            # count # of ben/harm (gene/cnv) per patient
            if not patient in patient_harm_ben_counts:
                patient_harm_ben_counts[patient] = OrderedDict()

                patient_harm_ben_counts[patient]['weighted'] = OrderedDict()
                patient_harm_ben_counts[patient]['unweighted'] = OrderedDict()

                patient_harm_ben_counts[patient]['weighted']['HARMFUL'] = 0
                patient_harm_ben_counts[patient]['weighted']['BENIGN'] = 0
                patient_harm_ben_counts[patient]['unweighted']['HARMFUL'] = 0
                patient_harm_ben_counts[patient]['unweighted']['BENIGN'] = 0
            patient_harm_ben_counts[patient]['weighted'][case_control] += weight
            patient_harm_ben_counts[patient]['unweighted'][case_control] += 1
        else:
            # HARDCODE: hardcode informative lines of arff file
            # remove this, was used before for testing
            if preline.count('{') > 1:
                arff_info_pre.append(['sample', line[2], 'cc', 'f', '0.0', 'none', 'phenotype', 'asdf'] + [preline.split('{')[2].split('}')[0]])
            else:
                arff_info_pre.append(['sample', line[2], 'cc', 'f', '0.0', 'none', 'phenotype', 'asdf'] + [1])

# pred
f = _gzip_open(args.pred_file, 'r')
# skip header
while 1:
    line = f.readline().strip().split()
    if len(line) == 5:
        break

# print class_values
out_pred_w_cnv_f = _gzip_open('out_pred_w_cnv%s.txt.gz' % args.out_file_postfix, 'w')
for line in f:
    line = line.strip().split()

    if len(line) < 4:
        continue
    
    if len(line) == 4:
        type_list = [int, str, str, float]
	ret = [b(a) for a,b in zip(line, type_list)]
        ret.insert(3, '-')
    elif len(line) == 5:
        type_list = [int, str, str, str, float]
	ret = [b(a) for a,b in zip(line, type_list)]

    (inst, actual, predicted, correct, conf) = ret

    new_correct = None
    if actual == '1:BENIGN':
        if arff_info_pre[inst][7] == 'max':
            new_correct = 'wrong'
        elif arff_info_pre[inst][7] == 'benign':
            new_correct = 'correct'
        elif arff_info_pre[inst][7] == 'notmax':
            new_correct = 'correct'
    elif actual == '2:HARMFUL':
        if arff_info_pre[inst][7] == 'max':
            new_correct = 'correct'
        elif arff_info_pre[inst][7] == 'benign':
            new_correct = 'wrong'
        elif arff_info_pre[inst][7] == 'notmax':
            new_correct = 'wrong'
    else:
        raise Error('not valid actual: ' % actual)
    
    # HARDCODE: index
    patient = arff_info_pre[inst][0]
    case_control = arff_info_pre[inst][2]
    arff_info_pre[inst] = ret + arff_info_pre[inst] + [new_correct, \
                                                       patient_harm_ben_counts[patient]['unweighted']['HARMFUL'], \
                                                       patient_harm_ben_counts[patient]['unweighted']['BENIGN'], \
                                                       patient_harm_ben_counts[patient]['weighted']['HARMFUL'], \
                                                       patient_harm_ben_counts[patient]['weighted']['BENIGN']]

    # debug
    print >> out_pred_w_cnv_f, '\t'.join(map(str, arff_info_pre[inst]))
out_pred_w_cnv_f.close()

new_actual_list = ['max', 'notmax', 'benign']

i = 1
arff_row_template_pre = copy.deepcopy(arff_info_pre[-1])
for act in new_actual_list:
    for pred in ['1:BENIGN', '2:HARMFUL']:
        arff_row_template = copy.deepcopy(arff_row_template_pre)
        arff_row_template[0] += i
        arff_row_template[12] = act
        arff_row_template[2] = pred
        arff_row_template[5] = 'Copy_%s_%s' % (act, pred)
        arff_row_template[6] = 'Copy_%s_%s' % (act, pred)
        if act == 'max':
            arff_row_template[7] = 'HARMFUL'
            arff_row_template[1] = '2:HARMFUL'
        else:
            arff_row_template[7] = 'BENIGN'
            arff_row_template[1] = '1:BENIGN'
        arff_info_pre.append(arff_row_template)
        i += 1

arff_info = DataFrame(arff_info_pre[1:], columns=['inst', 'actual', 'predicted', 'correct', 'conf', 'sample', 'cnv', 'case_control', 'dup', 'hposim', 'original_gene', 'hpo_term', 'infomax', 'maxgene', 'geneconf', 'genesim', 'weight', 'newcorrect', 'wnumharm', 'wnumben', 'uwnumharm', 'uwnumben'])
# print arff_info[-10:].to_string()

confusion_matrix = DataFrame(index=new_actual_list, columns=class_values)

for actual, predicted in itertools.product(new_actual_list, class_values):
    num = sum(arff_info[(arff_info.infomax == actual) & (arff_info.predicted == predicted)].weight)
    confusion_matrix[predicted].ix[actual] = num - 1

num_correct = confusion_matrix['2:HARMFUL'].ix['max'] + confusion_matrix['1:BENIGN'].ix['notmax'] + confusion_matrix['1:BENIGN'].ix['benign']
num_incorrect = confusion_matrix['1:BENIGN'].ix['max'] + confusion_matrix['2:HARMFUL'].ix['notmax'] + confusion_matrix['2:HARMFUL'].ix['benign']
total = num_correct + num_incorrect

print 'By genemax:'

print 'classified as -->\n|'
print confusion_matrix.to_string()

precision_harm_num = confusion_matrix['2:HARMFUL'].ix['max']
precision_harm_denom = confusion_matrix['2:HARMFUL'].ix['max'] + confusion_matrix['2:HARMFUL'].ix['notmax'] + confusion_matrix['2:HARMFUL'].ix['benign']
recall_harm_num = precision_harm_num
recall_harm_denom = confusion_matrix['2:HARMFUL'].ix['max'] + confusion_matrix['1:BENIGN'].ix['max']
precision_harm = precision_harm_num/float(precision_harm_denom)*100
recall_harm = recall_harm_num/float(recall_harm_denom)*100
fmeasure_harm_num = 2*precision_harm*recall_harm
fmeasure_harm_denom = precision_harm + recall_harm
fmeasure_harm = fmeasure_harm_num/fmeasure_harm_denom

print ' genemax Correctly Classified Instances:\t%s\t%s %%' % ("%s/%s" % (num_correct, total), num_correct/float(total)*100)
print ' genemax Incorrectly Classified Instances:\t%s\t%s %%' % ("%s/%s" % (num_incorrect, total), num_incorrect/float(total)*100)
print ' genemax Total Number of Instances:\t%s' % total
print ' genemax Precision on harmful:\t%s/%s\t%s %%' % (precision_harm_num, precision_harm_denom, precision_harm)
print ' genemax Recall on harmful:\t%s/%s\t%s %%' % (recall_harm_num, recall_harm_denom, recall_harm)
print ' genemax F-measure on harmful:\t%s/%s\t%s %%' % (fmeasure_harm_num, fmeasure_harm_denom, fmeasure_harm)




confusion_matrix = DataFrame(index=class_values, columns=class_values)

for actual, predicted in itertools.product(class_values, class_values):
    num = sum(arff_info[(arff_info.actual == actual) & (arff_info.predicted == predicted)].weight)
    confusion_matrix[predicted].ix[actual] = num - 1

# [predicted][actual]
confusion_matrix['2:HARMFUL'].ix['1:BENIGN'] -= 1
confusion_matrix['1:BENIGN'].ix['1:BENIGN'] -= 1

num_correct = confusion_matrix['2:HARMFUL'].ix['2:HARMFUL'] + confusion_matrix['1:BENIGN'].ix['1:BENIGN']
num_incorrect = confusion_matrix['1:BENIGN'].ix['2:HARMFUL'] + confusion_matrix['2:HARMFUL'].ix['1:BENIGN']
total = num_correct + num_incorrect

print '\n\nBy gene:'

print 'classified as -->\n|'
print confusion_matrix.to_string()

precision_harm_num = confusion_matrix['2:HARMFUL'].ix['2:HARMFUL']
precision_harm_denom = confusion_matrix['2:HARMFUL'].ix['1:BENIGN'] + confusion_matrix['2:HARMFUL'].ix['2:HARMFUL']
recall_harm_num = precision_harm_num
recall_harm_denom = confusion_matrix['2:HARMFUL'].ix['2:HARMFUL'] + confusion_matrix['1:BENIGN'].ix['2:HARMFUL']
precision_harm = precision_harm_num/float(precision_harm_denom)*100
recall_harm = recall_harm_num/float(recall_harm_denom)*100
fmeasure_harm_num = 2*precision_harm*recall_harm
fmeasure_harm_denom = precision_harm + recall_harm
fmeasure_harm = fmeasure_harm_num/fmeasure_harm_denom

print ' gene Correctly Classified Instances:\t%s\t%s %%' % ("%s/%s" % (num_correct, total), num_correct/float(total)*100)
print ' gene Incorrectly Classified Instances:\t%s\t%s %%' % ("%s/%s" % (num_incorrect, total), num_incorrect/float(total)*100)
print ' gene Total Number of Instances:\t%s' % total
print ' gene Precision on harmful:\t%s/%s\t%s %%' % (precision_harm_num, precision_harm_denom, precision_harm)
print ' gene Recall on harmful:\t%s/%s\t%s %%' % (recall_harm_num, recall_harm_denom, recall_harm)
print ' gene F-measure on harmful:\t%s/%s\t%s %%' % (fmeasure_harm_num, fmeasure_harm_denom, fmeasure_harm)




print '\n\nBy cnv:'
arff_pt = pivot_table(arff_info, values=[], cols=['actual', 'predicted'], rows=['cnv'], aggfunc=len, fill_value=0)

out_pt = _gzip_open('out_pt%s.txt.gz' % args.out_file_postfix, 'w')
print >> out_pt, arff_pt.to_string()
out_pt.close()


# print arff_pt[-10:]
# print arff_pt.to_string()
print 'Actual -> 2:HARMFUL\nPredicted -v'
harmful_pt = pivot_table((arff_pt > 0)['2:HARMFUL'], values=[], cols=[], rows=['1:BENIGN', '2:HARMFUL'], aggfunc=len)
# print harmful_pt
harmful_pt[0][0] -= 2
harmful_pt[0][1] -= 1
harmful_pt[1][0] -= 1
print harmful_pt

print

print 'Actual -> 1:BENIGN\nPredicted -v'
benign_pt = pivot_table((arff_pt > 0)['1:BENIGN'], values=[], cols=[], rows=['1:BENIGN', '2:HARMFUL'], aggfunc=len)
# print benign_pt
benign_pt[0][0] -= 2
benign_pt[0][1] -= 2
benign_pt[1][0] -= 2
print benign_pt

# print dict(harmful_pt)
num_correct = harmful_pt.get((False, True), 0) + harmful_pt.get((True, True), 0) + benign_pt.get((True, False), 0)
num_incorrect = harmful_pt.get((True, False), 0) + benign_pt.get((False, True), 0) + benign_pt.get((True, True), 0)
total = num_correct + num_incorrect

precision_harm_num = harmful_pt.get((False, True), 0) + harmful_pt.get((True, True), 0)
precision_harm_denom = harmful_pt.get((False, True), 0) + harmful_pt.get((True, True), 0) + benign_pt.get((False, True), 0) + benign_pt.get((True, True), 0)
recall_harm_num = precision_harm_num
recall_harm_denom = harmful_pt.get((False, True), 0) + harmful_pt.get((True, True), 0) + harmful_pt.get((True, False), 0)
precision_harm = precision_harm_num/float(precision_harm_denom)*100
recall_harm = recall_harm_num/float(recall_harm_denom)*100
fmeasure_harm_num = 2*precision_harm*recall_harm
fmeasure_harm_denom = precision_harm + recall_harm
fmeasure_harm = fmeasure_harm_num/fmeasure_harm_denom

print ' cnv Correctly Classified Instances:\t%s\t%s %%' % ("%s/%s" % (num_correct, total), num_correct/float(total)*100)
print ' cnv Incorrectly Classified Instances:\t%s\t%s %%' % ("%s/%s" % (num_incorrect, total), num_incorrect/float(total)*100)
print ' cnv Total Number of Instances:\t%s' % total
print ' cnv Precision on harmful:\t%s/%s\t%s %%' % (precision_harm_num, precision_harm_denom, precision_harm)
print ' cnv Recall on harmful:\t%s/%s\t%s %%' % (recall_harm_num, recall_harm_denom, recall_harm)
print ' cnv F-measure on harmful:\t%s/%s\t%s %%' % (fmeasure_harm_num, fmeasure_harm_denom, fmeasure_harm)
