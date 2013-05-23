#!/usr/bin/env python
#Author:  Alessandro L. Asoni <ale.luca.asoni@gmail.com>
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
from itertools import product
from libkmersvm import *

def get_kmer_list(k):
    kmerlist = []
    for x in product('ACGT', repeat=k):
        kmerlist.append(''.join(x))
    return kmerlist

def kpair2feat_map(kmers, homeopair, quiet):
    if not quiet:
        print "creating kmer-pair to feature id map ... "
    feature_id = 1
    kpair_id_dict = {}
    #variables to keep track of progress
    progress_count = 0
    total = len(kmers)
    for a in kmers:
        if not quiet and progress_count % math.ceil(total/100.0) == 0:
            p = (float(progress_count)/total)*100
            sys.stdout.write("\r%.2f%% progress. " % p )
            sys.stdout.flush()
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
        print
    return kpair_id_dict

def create_feature_vector( seq, features, k, dmin, dmax):
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


    return feature_vector

def spectrum_kernel( x1, x2 ):
    mag_x1 = math.sqrt(float(inner_product(x1,x1)))
    mag_x2 = math.sqrt(float(inner_product(x2,x2)))
    return inner_product(x1,x2)/(mag_x1*mag_x2)

def inner_product( x1, x2 ):
    in_product = 0
    for ( f_id, f_count) in x1.items():
        if f_id in x2:
            in_product += f_count * x2[f_id]

    return in_product

def create_kernel_matrix(n):
    M = numpy.zeros((n,n), dtype=float)
    return M

def update_kernel_matrix(M,n,m,d):
    M[n,m] = d
    return

def fill_kernel_matrix(seqs,M,features,k,dmin,dmax,quiet):
    if not quiet:
        print "filling kernel matrix ... "
    
    #memorize feature vecotrs
    feature_vectors = [{}] * len(seqs)

    #variables to keep trak of progress
    progress_count = 0
    total = len(seqs)
    for i in xrange(0,len(seqs)):
        if not quiet and progress_count % math.ceil(total/1000.0) == 0:
            p = (float(progress_count)/total)*100
            sys.stdout.write("\r%.2f%% progress. " % p)
            sys.stdout.flush()
        progress_count += 1
        for j in xrange(0,i+1):
            if not len(feature_vectors[i]) == 0:
                xi = feature_vectors[i]
            else:
                xi = create_feature_vector(seqs[i],features,k,dmin,dmax)
                feature_vectors[i] = xi

            if not len(feature_vectors[j]) == 0:
                xj = feature_vectors[j]
            else:
                xj = create_feature_vector(seqs[j],features,k,dmin,dmax)
                feature_vectors[j] = xj

            M[i,j] = spectrum_kernel( xi, xj )
            if i == j and M[i,j] == 0:
                print create_feature_vector(seqs[i],features,k,dmin,dmax)
                print create_feature_vector(seqs[j],features,k,dmin,dmax)
            
    if not quiet:
        print
    return M

def write_kernel_matrix(M,n,output):
    f = open(output, 'w')
    for i in xrange(0,n):
        s = " " + str(M[i,0])
        f.write(s)
        for j in xrange(1,i+1):
            s = " " + str(M[i,j])
            f.write(s)
        f.write("\n")
    return

def main(argv = sys.argv):
    usage = "Usage: %prog [options] POSITIVE_SEQ NEGATIVE_SEQ OUTPUT_NAME"
    desc  = "1. take two files(FASTA format) as input, 2. train an SVM and store the trained SVM weights"
    parser = optparse.OptionParser(usage=usage, description=desc)
    parser.add_option("-k", dest="kmerlen", type=int, help="set kmer length", default=6)
    parser.add_option("-d",  dest="dmin", type=int, help="set minimum distance between kmer pair", default=0)
    parser.add_option("-D",  dest="dmax", type=int, help="set maximum distance between kmer pair", default=50)
    parser.add_option("-H",  dest="homeopair", default=True, help="don't use duplicate kmer pair as feature (default=true)", action="store_false")
    parser.add_option("-q", dest="quiet", default=False, action="store_true", help="supress messages (default=false)")

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

    kpair_map = kpair2feat_map( get_kmer_list(options.kmerlen), options.homeopair, options.quiet )

    M = create_kernel_matrix(len(seqs))
    M = fill_kernel_matrix(seqs,M,kpair_map,options.kmerlen, options.dmin, options.dmax, options.quiet)
    write_kernel_matrix(M,len(seqs),outm)

if __name__ =='__main__': main()
