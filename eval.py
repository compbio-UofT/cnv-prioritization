#!/usr/bin/env python

import sys
import argparse
import os
import json
from collections import OrderedDict
import errno

def arg_parser():
    '''Parse the arguments of the script.
    '''
    parser = argparse.ArgumentParser(description='Run weka and evaluate the results.')
    parser.add_argument('other_params', help="Other parameters, see README.md.")
    parser.add_argument('--nordm', '-r', help="Set to not generate a new randomized data set.", action='store_true')
    parser.add_argument('--debug', '-d', help="Debug.", action='store_true')

    args = parser.parse_args()
    return args

def mkdir_p(path):
    '''http://stackoverflow.com/questions/600268/mkdir-p-functionality-in-python
    '''
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST:
            pass
        else: raise

def pbash_cd(path):
    print 'cd %s' % path
    os.chdir(path)

def get_params(file_name, file_index):
    '''Get variable parameters from the parameter file. (SGE array)
    '''
    with open(file_name) as params_file:
        header = params_file.readline().rstrip().split()
        for line_num, line in enumerate(params_file):
            if line_num == file_index:
                line = line.rstrip().split()
                return OrderedDict(zip(header, line))

def other_params(params):
    '''Get user params.
    '''
    return [p.split('=') for p in params.split(';')]

def param_tuple(params, keys):
    '''Extract a parameter tuple.
    '''
    return tuple([params[k] for k in keys])

def update_range_params(params, key):
    if key in params:
        if params[key] != 'None':
            start, end = [int(i) for i in params[key].split(',')]
            params['%s_start' % key] = start
            params['%s_end' % key] = end
        else:
            params[key] = []
    else:
        params[key] = []

def full_params(other):
    params = OrderedDict()
    extra_params = other_params(other)
    print extra_params
    params.update(extra_params)

    for key in ['gtrain', 'gtest', 'ctest']:
        update_range_params(params, key)
    
    if params['skip'] != 'false':
        print 'job skipped'
        sys.exit()
        
    print json.dumps(params, indent=4)

    return params

def bash(cmd):
    os.system(cmd)

def pbash(cmd):
    ppbash(cmd)
    bash(cmd)

def ppbash(cmd):
    print cmd
    sys.stdout.flush()

def set_classifier(params, classifier_type):
    params['classifier_type'] = params[classifier_type]

    if params['classifier_type'] == "dtree":
        params['weka_cmd'] = "%(java_prog)s -Xmx3g -cp /filer/tools/weka/weka-3-6-6/weka.jar:/filer/tools/libsvm/libsvm-3.11/java/libsvm.jar weka.classifiers.trees.J48" % params
        params['weka_params'] = "-C 0.25 -M 2"
    elif params['classifier_type'] == "svm":
        params['weka_cmd'] = "%(java_prog)s -Xmx3g -cp /filer/tools/weka/weka-3-6-6/weka.jar:/filer/tools/libsvm/libsvm-3.11/java/libsvm.jar weka.classifiers.functions.LibSVM" % params
        params['weka_params'] = "-K 2"
    elif params['classifier_type'] == "adaboost":
        params['weka_cmd'] = "%(java_prog)s -Xmx3g -cp /filer/tools/weka/weka-3-6-6/weka.jar:/filer/tools/libsvm/libsvm-3.11/java/libsvm.jar weka.classifiers.meta.AdaBoostM1" % params
        params['weka_params'] = "-W weka.classifiers.trees.DecisionStump"
    elif params['classifier_type'] == "rf":
        params['weka_cmd'] = "%(java_prog)s -Xmx3g -cp /filer/tools/weka/weka-3-6-6/weka.jar:/filer/tools/libsvm/libsvm-3.11/java/libsvm.jar weka.classifiers.trees.RandomForest" % params
        params['weka_params'] = ""
    elif params['classifier_type'] == "nb":
        params['weka_cmd'] = "%(java_prog)s -Xmx3g -cp /filer/tools/weka/weka-3-6-6/weka.jar:/filer/tools/libsvm/libsvm-3.11/java/libsvm.jar weka.classifiers.bayes.NaiveBayes" % params
        params['weka_params'] = ""

    return params

def set_files(params):
    # files
    params['constants'] = '/filer/jfoong/projects/cnvsk/data/raw/for_scripts/'
    params['scripts'] = '/filer/jfoong/home/jfoong/jfoong/projects/sandbox/cnv/'
    params['file_postfix'] = '%s_%s_%s' % param_tuple(params, ['function_number', 'top_x', 'weighted_gene_duplication'])
    params['file_postfix_more'] = '%s_%s_%s' % param_tuple(params, ['file_postfix', 'balance', 'array_id'])
    return params

def run_create_dirs_and_links(params, args):
    # create dirs/links
    mkdir_p(params['file_postfix_more'])
    pbash_cd(params['file_postfix_more'])
    for ftype in ['out', 'err']:
        try:
            os.symlink('../err/%s.%s' % (args.sge_task_id, ftype), '%s.%s' % (params['file_postfix_more'], ftype))
        except:
            pass
    print os.getcwd()

def run_debug_sge(params):
    # debug
    pbash('ulimit -a' % params)
    pbash('ls /proc | wc -l' % params)
    pbash('uname -a' % params)

def run_set_constants(params):
    
    # HARDCODE: TODO:
    params['go_start'] = 2
    params['go_end'] = 6
    params['go_end_p1'] = params['go_end'] + 1
    params['hp_start'] = 113
    params['hp_end'] = 135

    params['rmcmd'] = "weka_remove.py"

    # java_prog=/filer/tools-supa/java/6.11/jre1.6.0_11/bin/java
    params['java_prog'] = "/usr/bin/java"

    rmd_cmd_list = [
        "%(rmcmd)s -R %(go_start)s-%(go_end_p1)s" % params,
        "%(rmcmd)s -R 1,2,3" % params,
        "%(rmcmd)s -R 4" % params,
        "%(rmcmd)s -R 2" % params,
        "%(rmcmd)s -R 3" % params,
        "%(rmcmd)s -R 3,4" % params,
        "%(rmcmd)s -R 2,4" % params,
        "%(rmcmd)s -R 2,3" % params,
        ]

    arff_list_no_go = [
        'go_hp.arff.gz',
        'len_dgv_gene.arff.gz',
        'len_dgv.arff.gz',
        'dgv_gene.arff.gz',
        'len_gene.arff.gz',
        'len.arff.gz',
        'dgv.arff.gz',
        'gene.arff.gz',
    ]

    return params, rmd_cmd_list, arff_list_no_go

def run_into_iteration(iteration_count):
    iter_dir = "iteration_%s" % (iteration_count)
    mkdir_p(iter_dir)
    pbash_cd(iter_dir)


def run_create_arff_and_set_classifier(this_arff, this_rmd_cmd, params):
    # gene level
    if this_arff == 'go_hp.arff.gz':

        # create arff file for extra data
        pbash("(cat %(constants)s/weka_header.txt; ls ../%(folds)s*_%(iteration)s_cnvnbal.arff.gz | xargs gunzip -c | sort -k4,4 -k8,8) | gzip > test_fold%(folds)s_go_hp_extra.arff.gz" % params)
        params['crmcmd'] = this_rmd_cmd
        pbash("%(crmcmd)s -i test_fold%(folds)s_go_hp_extra.arff.gz -o test_fold%(folds)s_go_hp.arff.gz" % params)

        params = set_classifier(params, 'classifier_type_gene')

    # cnv level
    else:
        params = set_classifier(params, 'classifier_type_cnv')

    return params
        
def run_create_train_test_noharm(params):
    pbash('(ls job_fold*/cnv_len_dgv*_bal_* | head -n 1 | xargs zgrep -v ^{ ; zgrep -H ^{ job_fold*/cnv_len_dgv*_bal_* | cut -d : -f 2-) | gzip > train_noharm_len_dgv_gene.arff.gz')
    pbash('(ls *cnvfold*noharm*len_dgv_gene* | head -n 1 | xargs zgrep -v ^{; zgrep -H ^{ *cnvfold*noharm*len_dgv_gene* | cut -d : -f 2-) | gzip > test_noharm_len_dgv_gene.arff.gz')

    job_dir = 'job_noharm' % params
    mkdir_p(job_dir)
    pbash_cd(job_dir)

    # run basic weka, on nbal, stops here for curves
    pbash("%(weka_cmd)s %(weka_params)s -t ../train_noharm_len_dgv_gene.arff.gz -T ../test_noharm_len_dgv_gene.arff.gz -c first -d model.model | gzip > out.txt.gz" % params)
    pbash("%(weka_cmd)s -T ../test_noharm_len_dgv_gene.arff.gz -c first -l model.model -p 0 | gzip > out_pred.txt.gz" % params)
    pbash("calc_results_on_cnvs.py ../test_noharm_len_dgv_gene.arff.gz out_pred.txt.gz '' | gzip > calc_res.txt.gz" % params)
    pbash("noharm_merge.py %(constants)s/dbAll.gff out_pred_w_cnv.txt.gz | gzip > new_harm_in_harmless_patients.txt.gz" % params)

    pbash_cd('..')

def run_into_job(params):
    job_dir = 'job_%(fprefix)s_%(arff_file)s' % params
    mkdir_p(job_dir)
    pbash_cd(job_dir)

def set_main_file(this_arff, params, fold):
    params['fold'] = fold

    main_file = [
        "fold%(fold)s_go_hp_extra.arff.gz" % params,
        "cnvfold%(fold)s_len_dgv_gene.arff.gz" % params,
        ]

    params['arff_file'] = this_arff
    if this_arff == 'go_hp.arff.gz':
        print "gene"
        params['main_file_w_index'] = main_file[0]
        params['fprefix'] = 'fold%s' % params['fold']
    elif this_arff == 'len_dgv_gene.arff.gz':
        print "cnv, 1st"
        params['main_file_w_index'] = main_file[1]
        params['fprefix'] = 'cnvfold%s' % params['fold']
    else:
        print "cnv, rest"
        params['main_file_w_index'] = main_file[1]
        params['fprefix'] = 'cnvfold%s' % params['fold']

    return params

def run_create_train_test(this_arff, params):

    if this_arff == 'go_hp.arff.gz':
        # train, gene; bal
        pbash("(cat %(constants)s/weka_header.txt; ls ../*_%(iteration)s_cnvbal*gz | grep -v /%(folds)s_ | grep -v /%(fold)s_ | xargs gunzip -c | sort -k4,4 -k8,8) | gzip > train_bal_%(main_file_w_index)s" % params)
        pbash("(cat %(constants)s/weka_header.txt; ls ../*_%(iteration)s_cnvnbal*gz | grep -v /%(folds)s_ | grep -v /%(fold)s_ | xargs gunzip -c | sort -k4,4 -k8,8) | gzip > train_nbal_%(main_file_w_index)s" % params)
        # test, gene; bal and nbal
        pbash("(cat %(constants)s/weka_header.txt; ls ../*_%(iteration)s_cnvbal*gz | grep -v /%(folds)s_ | grep /%(fold)s_ | xargs gunzip -c | sort -k4,4 -k8,8) | gzip > test_bal_%(main_file_w_index)s" % params)
        pbash("(cat %(constants)s/weka_header.txt; ls ../*_%(iteration)s_cnvnbal*gz | grep -v /%(folds)s_ | grep /%(fold)s_ | xargs gunzip -c | sort -k4,4 -k8,8) | gzip > test_nbal_%(main_file_w_index)s" % params)
    elif this_arff == 'len_dgv_gene.arff.gz':
        # merge gene level annotations

        # train, cnv; bal and nbal
        pbash("(zgrep -v ^{ */cnv*dgv*_nbal_*%(fold)s.arff.gz; ls */*dgv* | grep -v %(fold)s.arff.gz | grep -v nbal | xargs zgrep -H ^{ | cut -d : -f 2- ) | gzip > train_bal_%(main_file_w_index)s" % params)
        pbash("(zgrep -v ^{ */cnv*dgv*%(fold)s.arff.gz; ls */*dgv* | grep -v %(fold)s.arff.gz | grep nbal | xargs zgrep -H ^{ | cut -d : -f 2- ) | gzip > train_nbal_%(main_file_w_index)s" % params)
        # test, cnv; bal and nbal
        pbash("ln -s */*dgv*_bal_*%(fold)s.arff.gz test_bal_%(main_file_w_index)s" % params)
        pbash("ln -s */*dgv*_nbal_*%(fold)s.arff.gz test_nbal_%(main_file_w_index)s" % params)

def run_remove_attributes(this_arff, this_rmd_cmd, params):

    # do not weka_remove if on the first arff type
    if this_arff != 'len_dgv_gene.arff.gz':
        for tt in ['train', 'test']:
            for bb in ['bal', 'nbal']:
                params['crmcmd'] = this_rmd_cmd
                params['tt'] = tt
                params['bb'] = bb
                # current
                pbash("%(crmcmd)s -i %(tt)s_%(bb)s_%(main_file_w_index)s -o %(tt)s_%(bb)s_%(fprefix)s_%(arff_file)s" % params)

                if this_arff != 'go_hp.arff.gz' and params['balance'] == 'bt_patient':
                    # remaining
                    pbash("%(crmcmd)s -i test_cnvfold%(fold)s_noharm_len_dgv_gene.arff.gz -o test_cnvfold%(fold)s_noharm_%(arff_file)s" % params)

    return params
                                    
def run_basic_weka(params):

    # run basic weka, on bal, goes on to cnv step
    pbash("%(weka_cmd)s %(weka_params)s -t ../train_bal_%(fprefix)s_%(arff_file)s -T ../test_bal_%(fprefix)s_%(arff_file)s -c first -d model.model | gzip > out.txt.gz" % params)
    pbash("%(weka_cmd)s -T ../test_bal_%(fprefix)s_%(arff_file)s -c first -l model.model -p 0 | gzip > out_pred.txt.gz" % params)
    pbash("calc_results_on_cnvs.py ../test_bal_%(fprefix)s_%(arff_file)s out_pred.txt.gz '' | gzip > calc_res.txt.gz" % params)

    # run basic weka, on nbal, stops here for curves
    pbash("%(weka_cmd)s -T ../test_nbal_%(fprefix)s_%(arff_file)s -c first -l model.model | gzip  > out_nbal.txt.gz" % params)
    pbash("%(weka_cmd)s -T ../test_nbal_%(fprefix)s_%(arff_file)s -c first -l model.model -p 0 | gzip  > out_pred_nbal.txt.gz" % params)
    pbash("calc_results_on_cnvs.py ../test_nbal_%(fprefix)s_%(arff_file)s out_pred_nbal.txt.gz _nbal | gzip > calc_res_nbal.txt.gz" % params)


# also run the gene model on the incorrect of the len_dgv.arff.gz model
def run_incorrect(this_arff, crmcmd, params):

    # create incorrect files
    if this_arff == 'len_dgv.arff.gz':
        pbash("dgv_arff.py %(constants)s/cnvs_w_dgv_overlap.txt.gz ../job_fold%(fold)s_go_hp.arff.gz/out_pt_nbal.txt.gz out_pred_w_cnv_nbal.txt.gz -i | gzip > incorrect_len_dgv_nbal_raw_%(fold)s.arff.gz" % params)

    # eval incorrect files
    if this_arff == 'gene.arff.gz':
        # pbash("(zgrep -v ^{ ../*/incorrect_len_dgv_nbal_raw_*%(fold)s.arff.gz; ls ../*/incorrect_len_dgv_nbal_raw_*.arff.gz | grep -v %(fold)s.arff.gz | xargs zgrep -H ^{ | cut -d : -f 2- ) | gzip > test_incorrect_len_dgv_nbal_%(fold)s_pre.arff.gz" % params)
        pbash('ln -s ../job_cnvfold%(fold)s_len_dgv.arff.gz/incorrect_len_dgv_nbal_raw_%(fold)s.arff.gz .' % params)

        # HARDCODE!
        params['crmcmd'] = crmcmd

        pbash("%(crmcmd)s -i incorrect_len_dgv_nbal_raw_%(fold)s.arff.gz -o test_incorrect_len_dgv_nbal_%(fold)s.arff.gz" % params)
        pbash("%(weka_cmd)s -T test_incorrect_len_dgv_nbal_%(fold)s.arff.gz -c first -l model.model | gzip  > out_incorrect.txt.gz" % params)
        pbash("%(weka_cmd)s -T test_incorrect_len_dgv_nbal_%(fold)s.arff.gz -c first -l model.model -p 0 | gzip  > out_pred_incorrect.txt.gz" % params)
        pbash("calc_results_on_cnvs.py test_incorrect_len_dgv_nbal_%(fold)s.arff.gz out_pred_incorrect.txt.gz _incorrect | gzip > calc_res_incorrect.txt.gz" % params)

    if params['classifier_type'] == "dtree":
        pbash("%(weka_cmd)s -T ../test_bal_%(fprefix)s_%(arff_file)s -c first -l model.model -g | gzip  > out_graph.txt.gz" % params)
        pbash("parse_j48graph.py out_graph.txt.gz %(constants)s/gene_ontology_ext.obo %(constants)s/human-phenotype-ontology.obo" % params)

def run_noharm(this_arff, params):
    # gene
    if this_arff == 'go_hp.arff.gz':
        # run weka on extra data (patients without a harmful cnv)
        if params['balance'] == 'bt_patient':
            pbash("%(weka_cmd)s -T ../test_fold%(folds)s_%(arff_file)s -c first -l model.model -p 0 | gzip > out_pred_noharm.txt.gz" % params)
            pbash("%(weka_cmd)s -T ../test_fold%(folds)s_%(arff_file)s -c first -l model.model | gzip > out_noharm.txt.gz" % params)
            pbash("calc_results_on_cnvs.py ../test_fold%(folds)s_%(arff_file)s out_pred_noharm.txt.gz _noharm | gzip > calc_res_noharm.txt.gz" % params)

        # if first arff type, create cnv arff
        pbash("dgv_arff.py %(constants)s/cnvs_w_dgv_overlap.txt.gz out_pt.txt.gz out_pred_w_cnv.txt.gz | gzip > cnv_len_dgv_gene_bal_%(fold)s.arff.gz" % params)
        pbash("dgv_arff.py %(constants)s/cnvs_w_dgv_overlap.txt.gz out_pt_nbal.txt.gz out_pred_w_cnv_nbal.txt.gz | gzip > cnv_len_dgv_gene_nbal_%(fold)s.arff.gz" % params)

        if params['balance'] == 'bt_patient':
            pbash("dgv_arff.py %(constants)s/cnvs_w_dgv_overlap.txt.gz out_pt_noharm.txt.gz out_pred_w_cnv_noharm.txt.gz | gzip > ../test_cnvfold%(fold)s_noharm_len_dgv_gene.arff.gz" % params)

    # cnv
    else:
        if params['balance'] == 'bt_patient':
            # run weka on extra data (patients without a harmful cnv)
            pbash("%(weka_cmd)s -T ../test_%(fprefix)s_noharm_%(arff_file)s -c first -l model.model -p 0 | gzip > out_pred_noharm.txt.gz" % params)
            pbash("calc_results_on_cnvs.py ../test_%(fprefix)s_noharm_%(arff_file)s out_pred_noharm.txt.gz _noharm | gzip > calc_res_noharm.txt.gz" % params)
            pbash("noharm_merge.py %(constants)s/dbAll.gff out_pred_w_cnv_noharm.txt.gz | gzip > new_harm_in_harmless_patients.txt.gz" % params)

def main():
    '''Main.
    '''

    # get python params
    args = arg_parser()
    # get sge params
    params = full_params(args.other_params)
    params = set_files(params)
    run_create_dirs_and_links(params, args)

    # generate randomized sets
    if not args.nordm:
        pbash('randomize_for_weka.py %(sim_file)s %(iteration)s %(file_postfix)s_%(balance)s %(constants)s/weka_header.txt %(balance)s %(folds)s %(remaining)s %(weighted_gene_duplication)s' % params)

    run_debug_sge(params)

    params, rmd_cmd_list, arff_list_no_go = run_set_constants(params)

    for iteration_count in range(1, int(params['iteration']) + 1):

        print 'iteration_count: %s' % iteration_count
        run_into_iteration(iteration_count)

        # rmd_cmd_list, 1 gene, 7 cnv
        for this_arff, this_rmd_cmd in zip(arff_list_no_go, rmd_cmd_list):

            print "this_arff: %s" % this_arff
            params = run_create_arff_and_set_classifier(this_arff, this_rmd_cmd, params)

            # only do gene except for blessed model
            if this_arff == 'go_hp.arff.gz' or \
                ((params['top_x'] == '4' or \
                   params['top_x'] == '10') and \
                 params['function_number'] == '2' and \
                 params['weighted_gene_duplication'] == 'sim'):

                # create overall train/test for noharm
                if this_arff == 'len_dgv_gene.arff.gz':
                    run_create_train_test_noharm(params)

                # cross fold
                if params['gtrain'] == []:

                    for fold in range(int(params['folds'])):
                        print "fold: %s" % fold

                        params = set_main_file(this_arff, params, fold)
                        run_create_train_test(this_arff, params)
                        params = run_remove_attributes(this_arff, this_rmd_cmd, params)
                        run_into_job(params)
                        run_basic_weka(params)
                        run_incorrect(this_arff, rmd_cmd_list[-1], params)
                        run_noharm(this_arff, params)
                        pbash_cd('..') # out of job

                # simple 3 sets
                else:
                    if this_arff == 'go_hp.arff.gz':

                        set_main_file(this_arff, params, 'gene')

                        # train, gene; bal and nbal
                        pbash("(cat %(constants)s/weka_header.txt; for gtrain in `seq %(gtrain_start)s %(gtrain_end)s`; do gunzip -c ../${gtrain}_*_%(iteration)s_cnvbal*gz; done | sort -k4,4 -k8,8) | gzip > train_bal_%(main_file_w_index)s" % params)
                        pbash("(cat %(constants)s/weka_header.txt; for gtrain in `seq %(gtrain_start)s %(gtrain_end)s`; do gunzip -c ../${gtrain}_*_%(iteration)s_cnvnbal*gz; done | sort -k4,4 -k8,8) | gzip > train_nbal_%(main_file_w_index)s" % params)
                    

                        # test, gene; bal and nbal
                        pbash("(cat %(constants)s/weka_header.txt; for gtest in `seq %(gtest_start)s %(gtest_end)s`; do gunzip -c ../${gtest}_*_%(iteration)s_cnvbal*gz; done | sort -k4,4 -k8,8) | gzip > test_bal_%(main_file_w_index)s" % params)
                        pbash("(cat %(constants)s/weka_header.txt; for gtest in `seq %(gtest_start)s %(gtest_end)s`; do gunzip -c ../${gtest}_*_%(iteration)s_cnvnbal*gz; done | sort -k4,4 -k8,8) | gzip > test_nbal_%(main_file_w_index)s" % params)

                        set_main_file(this_arff, params, 'cnv')

                        # train, cnv; bal and nbal
                        pbash("ln -s train_bal_foldgene_go_hp_extra.arff.gz train_bal_foldcnv_go_hp_extra.arff.gz" % params)
                        pbash("ln -s train_nbal_foldgene_go_hp_extra.arff.gz train_nbal_foldcnv_go_hp_extra.arff.gz " % params)
                    
                        # test, cnv; bal and nbal
                        pbash("(cat %(constants)s/weka_header.txt; for ctest in `seq %(ctest_start)s %(ctest_end)s`; do gunzip -c ../${ctest}_*_%(iteration)s_cnvbal*gz; done | sort -k4,4 -k8,8) | gzip > test_bal_%(main_file_w_index)s" % params)
                        pbash("(cat %(constants)s/weka_header.txt; for ctest in `seq %(ctest_start)s %(ctest_end)s`; do gunzip -c ../${ctest}_*_%(iteration)s_cnvnbal*gz; done | sort -k4,4 -k8,8) | gzip > test_nbal_%(main_file_w_index)s" % params)

                        run_remove_attributes(this_arff, this_rmd_cmd, params)

                        set_main_file(this_arff, params, 'gene')

                    set_main_file(this_arff, params, 'cnv')
                    run_remove_attributes(this_arff, this_rmd_cmd, params)
                    run_into_job(params)
                    run_basic_weka(params)

                    if this_arff == 'go_hp.arff.gz':
                        set_main_file(this_arff, params, 'cnv')
                        for bb in ['bal', 'nbal']:
                            for tt in ['train', 'test']:
                                params['tt'] = tt
                                params['bb'] = bb
                                pbash("%(weka_cmd)s -T ../%(tt)s_%(bb)s_foldcnv_go_hp.arff.gz -c first -l model.model | gzip  > out_c%(tt)sc%(bb)s.txt.gz" % params)
                                pbash("%(weka_cmd)s -T ../%(tt)s_%(bb)s_foldcnv_go_hp.arff.gz -c first -l model.model -p 0 | gzip  > out_pred_c%(tt)s%(bb)s.txt.gz" % params)
                                pbash("calc_results_on_cnvs.py ../%(tt)s_%(bb)s_foldcnv_go_hp.arff.gz out_pred_c%(tt)s%(bb)s.txt.gz _c%(tt)s%(bb)s | gzip > calc_res_c%(tt)s%(bb)s.txt.gz" % params)

                                pbash("dgv_arff.py %(constants)s/cnvs_w_dgv_overlap.txt.gz out_pt_c%(tt)s%(bb)s.txt.gz out_pred_w_cnv_c%(tt)s%(bb)s.txt.gz | gzip > ../%(tt)s_%(bb)s_cnvfoldcnv_len_dgv_gene.arff.gz" % params)

                    run_incorrect(this_arff, rmd_cmd_list[-1], params)
                    run_noharm(this_arff, params)

                    pbash_cd('..') # out of job

        pbash_cd('..') # out of iteration
                        
main()
