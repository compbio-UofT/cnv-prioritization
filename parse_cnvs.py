#!/usr/bin/env python

import csv
import pprint
import sys
import argparse
from collections import OrderedDict
import json
import copy

pp = pprint.PrettyPrinter(indent=0)

parser = argparse.ArgumentParser(description='Format input data into a gff file.')

parser.add_argument('csv_file_name', help="File for the original data (path + file).")
parser.add_argument('source', help="File label.")
parser.add_argument('case_list', help="List of strings to search for which indicate the cnv is clinically significiant.")
parser.add_argument('ctrl_list', help="List of strings to search for which indicate the cnv is a conrol cnv")
parser.add_argument('pheno_list', help="Phenotype list to filter on, '*' for all.")

args = parser.parse_args()

case_list = args.case_list.split('|')
ctrl_list = args.ctrl_list.split('|')
pheno_list = args.pheno_list.split('|')

# phenotypes to remove
OUT_LIST = ['IRRELEVANT', 'FAMILY_FOLLOW_UP', 'CHROMOSOME_ABNORMALIRY', 'CHROMOSOME_ABNORMALITY', \
            'OTHER',

            # added
            'Multiple Congenital Anomalies',
            'Congenital Anomaly',
            'Dysmorphic Features',
            'Parent of child with chromosome anomaly',
            'Congenital Heart Defect',
            'Family history of Autism',
            '? Prader-Willi Syndrome',
            'Epilepsy',
            'Reason Not Provided',
            'Family History',
            'Genomic array imbalance in family',

            # more
            'CONGENITAL_MALFORMATION',
            'MIM:176270',
            'MIM:213300',
            'MIM:608638',
            ]
RENAME_MAP = {
    'Developmental Delay' : 'HP:0001263',
    'Autism' : 'HP:0000717',
    'Facial Dysmorphism' : 'HP:0001999',
    'Cardiac Malformation' : 'HP:0002564',
    'Seizures' : 'HP:0001250',
    'Microcephaly' : 'HP:0000252',
    'Behavioural/Psychiatric Abnormality' : 'HP:0000708',
    'Cleft Palate' : 'HP:0000175',
    'Failure to Thrive' : 'HP:0001508',
    'Short Stature' : 'HP:0004322',
    'Intellectual Disability' : 'HP:0001249',
    'Learning Disability' : 'HP:0001328',
    'Ambiguous Genitalia' : 'HP:0000062',
    'Ventricular Septal Defect' : 'HP:0001629',
    'Hypospadias' : 'HP:0000047',
    'Stillbirth' : 'HP:0001624',
    'Increased risk of malignancy' : 'HP:0006741',
    'Amenorrhea' : 'HP:0000141',
    'Intrauterine growth restriction' : 'HP:0001511',
    'Speech delay' : 'HP:0002117',
    }

TRANSFORM_VAL_MAP = {
    'TYPE':{
    'Gain':'DUP',
    'Loss':'DEL'},
    }

def make_gff(variant, i, counts):
    '''Return a gff string.'''

    # True if the CNV is fine to print
    cnv_passes = True

    # Label case and controls
    if 'CLASSIFICATION' in variant:
        if variant['CLASSIFICATION'] in ctrl_list:
            cnvclass = 'CTL'
            phenotype_prefix = 'CONTROL_'
        elif variant['CLASSIFICATION'] in case_list:
            cnvclass = 'CNV'
            phenotype_prefix = ''
        else:
            # only print valid cases and controls
            cnv_passes = False
            cnvclass = 'ELSE'
            phenotype_prefix = 'ELSE_'

    # gff extra column
    extra = []
    for w in ['GENDER', 'SAMPLE_ID', 'TYPE', 'CYTOBAND', 'CLASSIFICATION']:
        if w in TRANSFORM_VAL_MAP and variant[w] in TRANSFORM_VAL_MAP[w]:
            to_append = TRANSFORM_VAL_MAP[w][variant[w]]
        else:
            to_append = variant[w]
        extra.append('%s=%s' % (w, to_append))

    # phenotypes extra column
    current_phenotype_list = []

    # only add 1 phenotype
    does_not_have_a_phenotype = True

    if type(variant['PHENOTYPE']) is list:
        phenotype_list = variant['PHENOTYPE']
    else:
        phenotype_list = [variant['PHENOTYPE']]

    for current_phenotype in phenotype_list:

        # TODO: use multiple phenotypes?
        if does_not_have_a_phenotype:

            # filter phenotypes by args
            if (pheno_list == ['*'] and (not current_phenotype in OUT_LIST)) or current_phenotype in pheno_list:

                if current_phenotype in RENAME_MAP:
                    current_phenotype = RENAME_MAP[current_phenotype]

                current_phenotype_list.append(phenotype_prefix + current_phenotype)

                if not current_phenotype in counts:
                    counts[current_phenotype] = OrderedDict()
                if not cnvclass in counts[current_phenotype]:
                    counts[current_phenotype][cnvclass] = 0

                counts[current_phenotype][cnvclass] += 1

                does_not_have_a_phenotype = False

    phenotype_str = "[" + ",".join(current_phenotype_list) + "]"            
            
    extra.append("PHENOTYPE=%s" % phenotype_str)

    # format main gff line
    line = ['chr' + variant['CHR'], args.source, "%s.%s.%08d" % (cnvclass, args.source, i), variant['START'], variant['STOP'], '.', '.', '.']
    # format last gff column
    line.append(";".join(extra))

    # only print the line if it includes the filtered phenotypes
    if current_phenotype_list == []:
        cnv_passes = False

    if cnv_passes:
    	print "\t".join(line)

    return counts

def change_header_120116(header):
    replace_dict = {
        'Date Scanned':'DATE_SCANNED',
        'Lab no.':'SAMPLE_ID',
        'Sex':'GENDER',
        'Clinical Indication':'PHENOTYPE',
        'Protocol':'PROTOCOL',
        'Genome Build':'GENOME_BUILD',
        'chromosome':'CHR',
        'start':'START',
        'stop':'STOP',
        'band':'CYTOBAND',
        'size':'SIZE',
        'max start':'MAX_START',
        'max stop':'MAX_STOP',
        'Max size':'MAX_SIZE',
        'No. probes':'NUM_PROBES',
        'Gain/Loss':'TYPE',
        'Classification':'CLASSIFICATION',
        'Inheritance':'INHERITANCE',
        'Log2 ratio':'LOG2_RATIO',
        'ISCN':'ISCN',
        }

    new_header = []
    for h in header:
        new_header.append(replace_dict.get(h, h))

    return new_header
    
firstline = file(args.csv_file_name).readline()
if len(firstline.split(",")) > len(firstline.split("\t")):
    csv_reader = csv.reader(open(args.csv_file_name, 'rb'), delimiter=',')
else:
    csv_reader = csv.reader(open(args.csv_file_name, 'rb'), delimiter='\t', quotechar='"')

i = -1

counts = OrderedDict()    

for row in csv_reader:
    if i == -1:
        header = row

        header = change_header_120116(header)

    else:
        variant = {}
        
        # transform line of variant info into a dict based on the header
        for j in range(len(row)):

            # handle ALL_GENES column marta said don't use this column, intersect with genes myself
            if header[j] == 'ALL_GENES':
                variant[header[j]] = [word.strip().split('#') for word in row[j][1:-1].split(",")]

            # handle list columns
            if len(row[j]) > 0 and row[j][0] == "[":
                variant[header[j]] = [word.strip() for word in row[j][1:-1].split(",")]

            # handle single value columns
            else:
                variant[header[j]] = row[j]

        # pp.pprint(variant)

        if not '180' in variant['SAMPLE_ID']:
            counts = make_gff(variant, i, counts)
        
    i += 1



count_summary = OrderedDict()
count_summary['VALID'] = 0
for phenotype, type_dict in counts.iteritems():
    for name, val in type_dict.iteritems():
        if not name in count_summary:
            count_summary[name] = 0
        count_summary[name] += counts[phenotype][name]

        if name in ['CNV', 'CTL']:
            count_summary['VALID'] += counts[phenotype][name]

f = open("%s_cnv_stats.txt" % args.source.split(".")[0], 'w')
print >> f, json.dumps(counts, indent=4)
print >> f, json.dumps(count_summary, indent=4)
f.close()

actual_all = {}
