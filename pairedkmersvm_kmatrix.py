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
from itertools import product
from libkmersvm import *


class Index(object):
    """ Substring index using hash table map """
    def __init__(self, t, ln=6):
        """ Create index, extracting substrings of length 'ln' """
        self.ln = ln
        self.index = {}
        for i in xrange(0, len(t)-ln+1):
            substr = t[i:i+ln]
            if substr not in self.index:
                self.index[substr] = []
            self.index[substr].append(i)

    def query(self, p):
        """ Return candidate alignments for p """
        return self.index.get(p[:self.ln]) or []

def kmerid2kmer(kmerid, kmerlen):
	"""convert integer kmerid to kmer string

	Arguments:
	kmerid -- integer, id of k-mer
	kmerlen -- integer, length of k-mer

	Return:
	kmer string
	"""

	nts = "ACGT"
	kmernts = []
	kmerid2 = kmerid

	for i in xrange(kmerlen):
		ntid = kmerid2 % 4
		kmernts.append(nts[ntid])
		kmerid2 = int((kmerid2-ntid)/4)

	return ''.join(reversed(kmernts))


def kmer2kmerid(kmer, kmerlen):
	"""convert kmer string to integer kmerid

	Arguments:
	kmerid -- integer, id of k-mer
	kmerlen -- integer, length of k-mer

	Return:
	id of k-mer
	"""

	nt2id = {'A':0, 'C':1, 'G':2, 'T':3}

	return reduce(lambda x, y: (4*x+y), [nt2id[x] for x in kmer])


def get_rcmap(kmerid, kmerlen):
	"""mapping kmerid to its reverse complement k-mer on-the-fly

	Arguments:
	kmerid -- integer, id of k-mer
	kmerlen -- integer, length of k-mer

	Return:
	integer kmerid after mapping to its reverse complement
	"""

	#1. get kmer from kmerid
	#2. get reverse complement kmer
	#3. get kmerid from revcomp kmer
	rckmerid = kmer2kmerid(revcomp(kmerid2kmer(kmerid, kmerlen)), kmerlen)

	if rckmerid < kmerid:
		return rckmerid

	return kmerid


def word_kernel( seq1, ind2, k, dmin, dmax ):
    """ Calculate kernel function applied to the inputs 
    Arguments:
        seq1 -- string, first sequence
        ind2 -- index, index for second sequences

    Return:
        kernel function of seq1 and seq2
    """
    dot_product = 0
    for i in xrange(0, len(seq1)-k+1):
        a = seq1[i:i+k]
        a_in_2 = ind2.query(a) + ind2.query(revcomp(a))
        if not a_in_2: continue

        left  = seq1[max(0, i-(dmax+k)):max(0, i-dmin)]
        right = seq1[min(len(seq1),i+k+dmin):min(len(seq1),i+dmax+k+k)]
        #checking for pairs to the left of 'a'
        for j in xrange(0, len(left)-k+1):
            b = left[j:j+k]
            b_in_2 = ind2.query(b) + ind2.query(revcomp(b))
            if not b_in_2: continue

            for m in a_in_2:
                for n in b_in_2:
                    if abs(m - n) <= dmax:
                        dot_product += 1

        #checking for pairs to the right of 'a'
        for j in xrange(0, len(right)-k+1):
            b = right[j:j+k]
            b_in_2 = ind2.query(b) + ind2.query(revcomp(b))
            if not b_in_2: continue

            for m in a_in_2:
                for n in b_in_2:
                    if abs(m - n) <= dmax:
                        dot_product += 1

    return dot_product


def create_kernel_matrix(n):
    """ Initialize kernel matrix 
    Arguments:
        n -- int, matrix dimension

    Return:
        the matrix 
    """
    M = numpy.zeros((n,n), dtype=float)
    return M


def fill_kernel_matrix(seqs,M,k,dmin,dmax,quiet):
    """ Fill kernel matrix 
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
        print "filling kernel matrix ... "
    
    #memorize index for each sequence
    indices = [None] * len(seqs)
    #memorize squared magnitudes
    square_mags = [0.0] * len(seqs)

    #variables to keep trak of progress
    progress_count = 0
    total = ( (len(seqs)+2)*(len(seqs)+1) )/2.0
    for i in xrange(0,len(seqs)):
        for j in xrange(i,-1,-1):
            if not quiet and progress_count % math.ceil(total/10000.0) == 0:
                p = (float(progress_count)/total)*100
                sys.stdout.write("\r%.2f%% progress. " % p)
                sys.stdout.flush()
            progress_count += 1
            
            if 

            if indices[j]:
                indj = indices[j]
                sqmj = square_mags[j]
            else:
                indj = Index(seqs[j],k)
                indeces[j] = indj

            krnl = word_kernel( seqs[i], indj, k, dmin, dmax)
            if i == j:
                M[i,j] = 1.0
            else 

            
            
    if not quiet:
        print

    return M


def write_kernel_matrix(M,n,output):
    """ write kernel matrix to file
    Arguments:
        M -- numpy matrix, kernel matrix
        n -- int, matrix dimension
        output -- string, name of output file
    Return:
        the matrix 
    """
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

    M = create_kernel_matrix(len(seqs))
    M = fill_kernel_matrix(seqs,M,options.kmerlen, options.dmin, options.dmax, options.quiet)
    write_kernel_matrix(M,len(seqs),outm)

if __name__ =='__main__': main()

