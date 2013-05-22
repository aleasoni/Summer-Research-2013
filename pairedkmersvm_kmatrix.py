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
from libkmersvm import *

def kpairs_in_seq(seq, k, dmin, dmax):
    """ Read a sequence. Return list of kmer pairs extracted 

    Arguments:
        seq -- string, DNA sequence
        k -- integer, kmer length
        dmin -- integer, mininmum distance b/w kmer pair
        dmax -- integer, maximum distance b/w kmer pair

    Return:
        list of kmer pairs (list of tuples)
        """
        kpairlist = []
        return kpairlist

def kpairs_in_window(seq, k, dmin, wsize, w):
    return

def create_kpair2feat_map(kpairlist):
    return


def create_feature_vector(label, seq, features, k, min_dist, max_dist):
    """ Create SVM-light format feature vector """
    return

def create_kernel_matrix(n,m):
    return

def update_kernel_matrix(D,n,m,d):
    return

def write_kernel_matrix(D):
    return

def main(argv = sys.argv):
    usage = "Usage: %prog [options] POSITIVE_SEQ NEGATIVE_SEQ"
    desc  = "1. take two files(FASTA format) as input, 2. train an SVM and store the trained SVM weights"
    parser = optparse.OptionParser(usage=usage, description=desc)
    parser.add_option("-k", dest="kmerlen", type=int, help="set kmer length", default=6)
    parser.add_option("-w", dest="wsize", type=int, help="set window size", default=100)
    parser.add_option("-d",  dest="dmin", type=int, help="set minimum distance between kmer pair", default=0)
    parser.add_option("-D",  dest="dmax", type=int, help="set maximum distance between kmer pair", default=0)

    (options, args) = parser.parse_args()

    if len(args) == 0:
        parser.print_help()
        sys.exit(0)

    if len(args) != 2:
        parser.error("incorrect number of arguments")
        parser.print_help()
        sys.exit(0)

    seqs_pos, sids_pos = read_fastafile(posf)
    seqs_neg, sids_neg = read_fastafile(negf)
    npos = len(seqs_pos)
    nneg = len(seqs_neg)
    seqs = seqs_pos + seqs_neg
    sids = sids_pos + sids_neg

if __name__ =='__main__': main()

