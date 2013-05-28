#!/usr/bin/env python
"""
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import sys
import string
import os
import optparse
import numpy
import math
import random
from subprocess import call
from itertools import product
from libkmersvm import *


def get_kmer_list(k):
    """ Genrate list of all possible kmers of length k 
    Arguments:
        k -- int, kmer length

    Return:
        list of all possible kmers of length k
    """
    kmerlist = []
    for x in product('ACGT', repeat=k):
        kmerlist.append(''.join(x))
    return kmerlist


def kpair2feat_map(kmers, homeopair, quiet):
    """ Create kmair-pair to feature id map 
    Arguments:
        kmers -- list of strings, list of all single kmers
        homeopair -- boolean, if true a pair of equal kmers is used as feature
        quiet -- boolean, if true it suppresses messages

    Return:
        dictionary of kmer-pairs (a,b) to feature id number
    """
    if not quiet:
        sys.stderr.write("creating kmer-pair to feature id map ...\r")
        sys.stderr.write("\n")
    feature_id = 1
    kpair_id_dict = {}
    #variables to keep track of progress
    progress_count = 0
    total = len(kmers)
    for a in kmers:
        if not quiet and progress_count % math.ceil(total/100.0) == 0:
            p = (float(progress_count)/total)*100
            sys.stderr.write("\r%.2f%% progress. " % p )
            sys.stderr.flush()
        progress_count += 1

        for b in kmers:
            if not homeopair:
                if a == revcomp(b): continue
                if a == b:          continue

            if (a,b) in kpair_id_dict: continue
            if (b,a) in kpair_id_dict: continue

            kpair_id_dict[ (a,b) ] = feature_id
            kpair_id_dict[ (a,revcomp(b)) ] = feature_id
            kpair_id_dict[ (revcomp(a),b) ] = feature_id
            kpair_id_dict[ (revcomp(a),revcomp(b)) ] = feature_id

            feature_id += 1
    if not quiet:
        sys.stderr.write("\n")
    return kpair_id_dict


def create_feature_vector( seq, pn, features, k, dmin, dmax):
    """ Create SVM-light format feature vector 
    Arguments:
        seq -- string, DNA sequence
        features -- dictionary, kmer pair to feature id map
        k -- integer, kmer length
        dmin -- integer, mininmum distance b/w kmer pair
        dmax -- integer, maximum distance b/w kmer pair

    Return:
        feature vector for seq
    """
    feature_vector = {}
    for i in xrange(0, len(seq)-k+1):
        a = seq[i:i+k]
        left  = seq[max(0, i-(dmax+k)):max(0, i-dmin)]
        right = seq[min(len(seq),i+k+dmin):min(len(seq),i+dmax+k+k)]
        #checking for pairs to the left of 'a'
        for j in xrange(0, len(left)-k+1):
            b = left[j:j+k]
            if (a,b) in features:
                feature_id = features[(a,b)]
            elif (b,a) in features:
                feature_id = features[(b,a)]
            else: continue

            if feature_id not in feature_vector:
                feature_vector[feature_id] = 0
            feature_vector[feature_id] += 1

        #checking for pairs to the right of 'a'
        for j in xrange(0, len(right)-k+1):
            b = right[j:j+k]
            if (a,b) in features:
                feature_id = features[(a,b)]
            elif (b,a) in features:
                feature_id = features[(b,a)]
            else: continue

            if feature_id not in feature_vector:
                feature_vector[feature_id] = 0
            feature_vector[feature_id] += 1

    feature_list = []
    for (feature, count) in feature_vector.items():
        feature_list.append((feature,count))
    feature_list.sort(key=lambda tup: tup[0])
    vector = ''.join(str(pn))
    for (feature,count) in feature_list:
        vector += " " + str(feature) + ":" + str(count)

    return vector


def write_feature_vectors(seqs,labs,k,dmin,dmax,quiet,output):
    """ write feature vectors
    Arguments:
        seqs -- list of strings, DNA sequences 
        M -- numpy matrix, kernel matrix
        k -- integer, kmer length
        dmin -- integer, mininmum distance b/w kmer pair
        dmax -- integer, maximum distance b/w kmer pair
        quiet -- boolean, if true it suppresses messages

    Return:
        the matrix 
    """
    if not quiet:
        sys.stderr.write("writing feature vectors to file ...\r")
        sys.stderr.write("\n")

    f = open(output, 'w')
    
    #variables to keep trak of progress
    progress_count = 0
    total = len(seqs)
    for i in xrange(0, total):
        if not quiet and progress_count % math.ceil(total/10000.0) == 0:
            p = (float(progress_count)/total)*100
            sys.stderr.write("\r%.2f%% progress. " % p)
            sys.stderr.flush()
        progress_count += 1

        f.write(create_feature_vector( seqs[i], labs[i], kpair_map, k, dmin, dmax ))
        f.write("\n")

    f.close()

    if not quiet:
        sys.stderr.write("\n")

    return 


def save_predictions(output, preds, cvs):
    """save prediction 
    """
    f = open(output, 'w')
    f.write('\t'.join(["#seq_id", "SVM score", "label", "NCV"]) + "\n")
    for i in xrange(len(preds)): 
        f.write('\t'.join([preds[i][1], str(preds[i][2]), str(preds[i][3]), str(cvs[i])]) + "\n")
    f.close()


def generate_cv_list(ncv, n1, n2):
    """generate the N-fold cross validation list

    Arguments:
    ncv -- integer, number of cross-validation
    n1 -- integer, number of positives
    n2 -- integer, number of negatives

    Return:
    a list of N-fold cross validation
    """

    shuffled_idx_list1 = range(n1)
    shuffled_idx_list2 = range(n1,n1+n2)

    random.shuffle(shuffled_idx_list1)
    random.shuffle(shuffled_idx_list2)

    shuffled_idx_list = shuffled_idx_list1 + shuffled_idx_list2

    idx = 0
    icv = 0
    cv = [0] * (n1+n2)
    while(idx < (n1+n2)):
        cv[shuffled_idx_list[idx]] = icv

        idx += 1
        icv += 1
        if icv == ncv:
            icv = 0

    return cv


def split_cv_list(cvlist, icv, data):
    """split data into training and test based on cross-validation list

    Arguments:
    cvlist -- list, cross-validation list
    icv -- integer, corss-validation set of interest
    data -- list, data set to be splitted

    Return:
    a list of training set and a list of test set
    """

    tr_data = []
    te_data = []

    for i in xrange(len(data)):
        if cvlist[i] == icv:
            te_data.append(data[i])
        else:
            tr_data.append(data[i])
    
    return tr_data, te_data

def svm_learn(seqs, labs, icv, options):
    cv_train = "cv"+str(icv)+".train"
    cv_model = "cv"+str(icv)+".model"
    write_feature_vectors( seqs, labs, options.kmerlen, options.dmin, options.dmax, options.quiet, cv_train )
    cmd = "svm_learn " + cv_train + " " + cv_model
    
    if not options.quiet:
        sys.stderr.write("executing: " + cmd + " ...\n")

    with open(os.devnull, "w") as fnull:
        result = call(["svm_learn", cv_train, cv_model], stdout=fnull, stderr=fnull)

    cmd = "rm " + cv_train
    if not options.quiet:
        sys.stderr.write("executing: " + cmd + "\n")

    call(["rm", cv_train])

    return cv_model


def svm_classify(seqs_te, labs_te, icv, svm_cv, options):
    cv_pred = []
    cv_test = "cv"+str(icv)+".test"
    cv_model = svm_cv
    cv_output = "cv"+str(icv)+".output"
    write_feature_vectors( seqs_te, labs_te, options.kmerlen, options.dmin, options.dmax, options.quiet, cv_test )
    cmd = "svm_classify " + cv_test + " " + cv_model + " " + cv_output
    if not options.quiet:
        sys.stderr.write("executing: " + cmd + "\n")

    with open(os.devnull, "w") as fnull:
        result = call(["svm_classify", cv_test, cv_model, cv_output], stdout=fnull, stderr=fnull)

    f = open(cv_output, 'r')
    for line in f:
        cv_pred.append(float(line))

    cmd = "rm " + cv_output + " " + cv_test + " " + cv_model
    if not options.quiet:
        sys.stderr.write("executing: " + cmd + "\n")

    call(["rm", cv_output, cv_test, cv_model])

    return cv_pred


def main(argv = sys.argv):
    usage = "Usage: %prog [options] POSITIVE_SEQ NEGATIVE_SEQ OUTPUT_NAME"
    desc  = "1. take two files(FASTA format) as input, \
             2. train an SVM and store the trained SVM weights"
    parser = optparse.OptionParser(usage=usage, description=desc)
    parser.add_option("-k", dest="kmerlen", type=int, \
            help="set kmer length", default=6)
    parser.add_option("-d", dest="dmin", type=int, \
            help="set minimum distance between kmer pair", default=0)
    parser.add_option("-D", dest="dmax", type=int, \
            help="set maximum distance between kmer pair", default=50)
    parser.add_option("-H", dest="homeopair", default=True, \
            help="don't use duplicate kmer pair as feature", action="store_false")
    parser.add_option("-q", dest="quiet", default=False, action="store_true", \
            help="supress messages (default=false)")
    parser.add_option("-v", dest="ncv", type="int", default=0, \
            help="if set, it will perform N-fold cross-validation and generate a prediction file (default = 0)")
    (options, args) = parser.parse_args()

    if len(args) == 0:
        parser.print_help()
        sys.exit(0)

    if len(args) != 3:
        parser.error("incorrect number of arguments")
        parser.print_help()
        sys.exit(0)

    posf = args[0]
    negf = args[1]
    outm = args[2]

    seqs_pos, sids_pos = read_fastafile(posf)
    seqs_neg, sids_neg = read_fastafile(negf)
    npos = len(seqs_pos)
    nneg = len(seqs_neg)
    seqs = seqs_pos + seqs_neg
    sids = sids_pos + sids_neg

    #generate labels
    labels = [1]*npos + [-1]*nneg

    #global variable
    global kpair_map
    kpair_map = kpair2feat_map( get_kmer_list(options.kmerlen), options.homeopair, options.quiet )

    if options.ncv > 0:
        if options.quiet == False:
            sys.stderr.write('cross-validation...\n')

        cvlist = generate_cv_list(options.ncv, npos, nneg)
        labels_cv = []
        preds_cv = []
        sids_cv = []
        indices_cv = []
        for icv in xrange(options.ncv):
            if options.quiet == False:
                sys.stderr.write("running cross validation number " + str(icv+1)  + " ...\n")
            #split data into training and test set
            seqs_tr, seqs_te = split_cv_list(cvlist, icv, seqs) 
            labs_tr, labs_te = split_cv_list(cvlist, icv, labels)
            sids_tr, sids_te = split_cv_list(cvlist, icv, sids)
            indices_tr, indices_te = split_cv_list(cvlist, icv, range(len(seqs)))

            #train SVM
            svm_cv = svm_learn(seqs_tr, labs_tr, icv, options)

            #test on current cv round
            preds_cv = preds_cv + svm_classify(seqs_te, labs_te, icv, svm_cv, options)
            
            labels_cv = labels_cv + labs_te
            sids_cv = sids_cv + sids_te
            indices_cv = indices_cv + indices_te

        output_cvpred = outm + "_cvpred.out"
        prediction_results = sorted(zip(indices_cv, sids_cv, preds_cv, labels_cv), key=lambda p: p[0])
        save_predictions(output_cvpred, prediction_results, cvlist)

    else:
        write_feature_vectors(seqs,labels,options.kmerlen,options.dmin,options.dmax,options.quiet,outm)

if __name__ =='__main__': main()

