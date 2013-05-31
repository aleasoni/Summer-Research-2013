#!/usr/bin/python
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
import re
from subprocess import call, Popen, PIPE
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

def kmer2feat_map(kmers):
    """ Genrate single kmers to feature id map
    Arguments:
    kmers -- list, all kmers

    Return:
    dictionary of kmer to feature id
    """
    feature_id = 1

    kmer_id_dict = {}
    for kmer in kmers:
        if kmer in kmer_id_dict: continue

        kmer_id_dict[kmer] = feature_id
        kmer_id_dict[revcomp(kmer)] = feature_id

        feature_id += 1

    return kmer_id_dict


def kpair2feat_map(kmers, homeopair, counts, direction, quiet):
    """ Create kmair-pair to feature id map 
    Arguments:
    kmers -- list of strings, list of all single kmers
    homeopair -- boolean, if true a pair of equal kmers is used as feature
    counts -- boolean, if true use single kmer counts as features
    quiet -- boolean, if true it suppresses messages

    Return:
    dictionary of kmer-pairs (a,b) to feature id number
    """
    if not quiet:
        sys.stderr.write("Creating kmer-pair to feature id map...\n")

    feature_id = 1
    if counts:
        feature_id += max(kmer_map.itervalues())

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
            kpair_id_dict[ (revcomp(a),revcomp(b)) ] = feature_id

            if direction:
                feature_id += 1
                kpair_id_dict[ (a,revcomp(b)) ] = feature_id
                kpair_id_dict[ (revcomp(a),b) ] = feature_id
            else:
                kpair_id_dict[ (a,revcomp(b)) ] = feature_id
                kpair_id_dict[ (revcomp(a),b) ] = feature_id

            feature_id += 1
    if not quiet:
        sys.stderr.write("\n")
    return kpair_id_dict


def create_feature_vector( seq, pn, k, dmin, dmax, spectrum, counts):
    """ Create SVM-light format feature vector 
    Arguments:
    seq -- string, DNA sequence
    pn -- integer, pos/neg - [1/-1]
    k -- integer, kmer length
    dmin -- integer, mininmum distance b/w kmer pair
    dmax -- integer, maximum distance b/w kmer pair
    spectrum -- boolean, return normalized vector
    counts -- boolean, use single kmer counts as features

    Return:
    feature vector for seq
    """
    feature_vector = {}

    if counts:
        for i in xrange(0, len(seq)-k+1):
            kmer = seq[i:i+k]
            feature_id = kmer_map[kmer]

            if feature_id not in feature_vector:
                feature_vector[feature_id] = 0
            feature_vector[feature_id] += 1

    for i in xrange(0, len(seq)-k-k-dmax+1):
        a = seq[i:i+k]
        window = seq[min(len(seq),i+k+dmin):min(len(seq),i+dmax+k+k)]
        #checking for pairs to the right of 'a'
        for j in xrange(0, len(window)-k+1):
            b = window[j:j+k]
            if (a,b) in kpair_map:
                feature_id = kpair_map[(a,b)]
            elif (b,a) in kpair_map:
                feature_id = kpair_map[(b,a)]
            else: continue

            if feature_id not in feature_vector:
                feature_vector[feature_id] = 0
            feature_vector[feature_id] += 1

    feature_list = []
    mag = 0.0
    for (feature, count) in feature_vector.items():
        feature_list.append((feature,count))
        if spectrum:
            mag += count**2

    mag = math.sqrt(mag)
    feature_list.sort(key=lambda tup: tup[0])
    vector = ''.join(str(pn))
    for (feature,count) in feature_list:
        if spectrum:
            vector += " " + str(feature) + ":" + str(float(count)/mag)
        else: 
            vector += " " + str(feature) + ":" + str(count)

    return vector


def write_feature_vectors(seqs,labs,k,dmin,dmax,quiet,output,spectrum,counts):
    """ write feature vectors
    Arguments:
    seqs -- list of strings, DNA sequences 
    labs -- list, labels
    k -- integer, kmer length
    dmin -- integer, mininmum distance b/w kmer pair
    dmax -- integer, maximum distance b/w kmer pair
    quiet -- boolean, if true it suppresses messages
    output -- string, name of output file
    spectrum -- boolean, normalize feature vector
    counts -- boolean, use single kmer counts

    """
    if not quiet:
        sys.stderr.write("Writing feature vectors to file...\n")

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

        f.write(create_feature_vector( seqs[i], labs[i], k, dmin, dmax, spectrum, counts ))
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
    """train svm by calling svm_learn (from svm_light package)

    Arguments:
    seqs -- list, sequences 
    labs -- list, labels
    icv -- integer, corss-validation set of interest
    options -- object containing option data 

    Return:
    name of model file for the trained svm
    """
    cv_train = "cv"+str(icv)+ "_" + str(options.dmin) + "_" + str(options.dmax) + ".train"
    cv_model = "cv"+str(icv)+ "_" + str(options.dmin) + "_" + str(options.dmax) + ".model"
    write_feature_vectors( seqs, labs, options.kmerlen, options.dmin, options.dmax,\
                           options.quiet, cv_train, options.spectrum, options.counts )
    cmd = "svm_learn " + cv_train + " " + cv_model
    
    if not options.quiet:
        sys.stderr.write("Executing: " + cmd + "...\n")

    proc = Popen(["svm_learn", cv_train, cv_model], stdout=PIPE)

    if not options.quiet:
        p1 = re.compile("Reading*")
        p2 = re.compile("Setting default*")
        p3 = re.compile("Optimization finished*")
        p4 = re.compile("Number of*")
        while True:
            line = proc.stdout.readline()
            if line != '':
                m = p1.match(line)
                if m:
                    sys.stderr.write("[svm_learn]: done scanning examples.\n")
                m = p2.match(line)
                if m:
                    sys.stderr.write("[svm_learn]: done reading examples into memory.\n")
                m = p3.match(line)
                if m:
                    sys.stderr.write("[svm_learn]: done optimizing...")

                m = p4.match(line)
                if m:
                    sys.stderr.write( line.rstrip() + ".\n")
            else:
                break

    cmd = "rm " + cv_train
    if not options.quiet:
        sys.stderr.write("Executing: " + cmd + "\n")

    call(["rm", cv_train])

    return cv_model


def svm_classify(seqs_te, labs_te, icv, svm_cv, options):
    """test svm by calling svm_classify (from svm_light package)

    Arguments:
    seqs_te -- list, test sequences 
    labs_te -- list, test labels
    icv -- integer, corss-validation set of interest
    svm_cv -- string, name of trained svm model file
    options -- object containing option data 

    Return:
    name of model file for the trained svm
    """
    cv_pred = []
    cv_test = "cv"+str(icv)+ "_" + str(options.dmin) + "_" + str(options.dmax) + ".test"
    cv_model = svm_cv
    cv_output = "cv"+str(icv)+ "_" + str(options.dmin) + "_" + str(options.dmax) + ".output"
    write_feature_vectors( seqs_te, labs_te, options.kmerlen, options.dmin, options.dmax,\
                           options.quiet, cv_test, options.spectrum, options.counts )
    cmd = "svm_classify " + cv_test + " " + cv_model + " " + cv_output
    if not options.quiet:
        sys.stderr.write("Executing: " + cmd + "...\n")

    proc = Popen(["svm_classify", cv_test, cv_model, cv_output], stdout=PIPE)

    if not options.quiet:
        p1 = re.compile("Classifying test*")
        p2 = re.compile("Runtime \(*")
        p3 = re.compile("Accuracy \(*")
        while True:
            line = proc.stdout.readline()
            if line != '':
                m = p1.match(line)
                if m:
                    sys.stderr.write("[svm_classify]: done reading model.\n")

                m = p2.match(line)
                if m:
                    sys.stderr.write("[svm_classify]: done classifying examples...")

                m = p3.match(line)
                if m:
                    sys.stderr.write(line.rstrip() + ".\n")
            else:
                break

    f = open(cv_output, 'r')
    for line in f:
        cv_pred.append(float(line))

    cmd = "rm " + cv_output + " " + cv_test + " " + cv_model
    if not options.quiet:
        sys.stderr.write("Executing: " + cmd + "\n")

    call(["rm", cv_output, cv_test, cv_model])

    return cv_pred


def main(argv = sys.argv):
    usage = "Usage: %prog [options] POSITIVE_SEQ NEGATIVE_SEQ OUTPUT_NAME"
    parser = optparse.OptionParser(usage=usage)
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

    parser.add_option("-c", dest="counts", default=False, \
            help="use single kmer counts as features (default=false)", action="store_true")

    parser.add_option("-s", dest="spectrum", default=False, \
			help="set the type of kernel to spectrum (default=false)", action="store_true")

    parser.add_option("-x", dest="direct", default=False, \
			help="use distinct features for distinct revcomp directions (default=false)", action="store_true")

    (options, args) = parser.parse_args()


    if len(args) == 0:
        parser.print_help()
        sys.exit(0)

    if len(args) != 3:
        parser.error("incorrect number of arguments")
        parser.print_help()
        sys.exit(0)

    if options.dmax < options.dmin + options.kmerlen:
        sys.stderr.write("error: dmax must be >= (dmin + kmerlen)\n")
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

    #list of kmers of length kmerlen
    kmers = get_kmer_list(options.kmerlen)

    if options.counts:
        #global variable
        global kmer_map
        kmer_map = kmer2feat_map(kmers)

    #global variable
    global kpair_map
    kpair_map = kpair2feat_map( kmers, options.homeopair, options.counts, options.direct, options.quiet )

    if options.ncv > 0:

        cvlist = generate_cv_list(options.ncv, npos, nneg)
        labels_cv = []
        preds_cv = []
        sids_cv = []
        indices_cv = []
        for icv in xrange(options.ncv):
            if options.quiet == False:
                sys.stderr.write("--- CROSS VALIDATION  " + str(icv+1)  + " ---\n")
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
        write_feature_vectors(seqs,labels,options.kmerlen,options.dmin,options.dmax,options.quiet,outm, options.spectrum, options.counts)

if __name__ =='__main__': main()

