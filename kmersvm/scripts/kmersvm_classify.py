#!/usr/bin/python
"""
	kmersvm_classify.py; classify sequences using SVM
	Copyright (C) 2011 Dongwon Lee

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
import numpy
import optparse

from libkmersvm import *

"""
global variables
"""
g_kmer2id = {}


class Parameters:
	def __init__(self, kernel=None, kmerlen=None, kmerlen2=None, bias=None, A=None, B=None):
		self.kernel = kernel 
		self.kmerlen = kmerlen
		self.kmerlen2 = kmerlen2 
		self.bias = bias
		self.A = A
		self.B = B


def read_svmwfile_wsk(filename):
	"""read SVM weight file generated by kmersvm_train.py

	Arguments:
	filename -- string, name of the SVM weight file

	Return:
	list of SVM weights
	an object of Parameters class 
	"""

	try:
		f = open(filename, 'r')
		lines = f.readlines()
		f.close()

	except IOError, (errno, strerror):
		print "I/O error(%d): %s" % (errno, strerror)
		sys.exit(0)

	kmer_svmw_dict = {}
	params = Parameters()

	for line in lines:
		#header lines
		if line[0] == '#':
			#if this line contains '=', that should be evaluated as a parameter
			if line.find('=') > 0:
				name, value = line[1:].split('=')
				vars(params)[name] = value
		else:
			s = line.split()
			kmerlen = len(s[0])
			if kmerlen not in kmer_svmw_dict:
				kmer_svmw_dict[kmerlen] = {}

			kmer_svmw_dict[kmerlen][s[0]] = float(s[2])

	#type casting of parameters
	params.kernel = int(params.kernel)
	params.kmerlen = int(params.kmerlen)
	if params.kernel == 1:
		params.kmerlen2 = params.kmerlen
	else:
		params.kmerlen2 = int(params.kmerlen2)
	params.bias = float(params.bias)
	params.A = float(params.A)
	params.B = float(params.B)

	#set global variable
	global g_kmer2id
	for k in range(params.kmerlen, params.kmerlen2+1):
		kmers = generate_kmers(k)
		rcmap = generate_rcmap_table(k, kmers)
		for i in xrange(len(kmers)): 
			g_kmer2id[kmers[i]] = rcmap[i]
	
	#create numpy arrays of svm weights
	svmw_list = []
	for k in range(params.kmerlen, params.kmerlen2+1):
		svmw = [0]*(2**(2*k))

		for kmer in kmer_svmw_dict[k].keys():
			svmw[g_kmer2id[kmer]] = kmer_svmw_dict[k][kmer]

		svmw_list.append(numpy.array(svmw, numpy.double))

	return svmw_list, params


def score_seq(s, svmw, kmerlen):
	"""calculate SVM score of given sequence using single set of svm weights

	Arguments:
	s -- string, DNA sequence
	svmw -- numpy array, SVM weights 
	kmerlen -- integer, length of k-mer of SVM weight

	Return:
	SVM score
	"""
	kmer2id = g_kmer2id
	x = [0]*(2**(2*kmerlen))
	for j in xrange(len(s)-kmerlen+1):
		x[ kmer2id[s[j:j+kmerlen]] ] += 1

	x = numpy.array(x, numpy.double)
	score_norm = numpy.dot(svmw, x)/numpy.sqrt(numpy.sum(x**2))

	return score_norm


def score_seq_wsk(s, svmwlist, kmerlen_start, kmerlen_end):
	"""calculate svm score of given sequence with multiple sets of svm weights

	Arguments:
	svmwlist -- list, SVM weights
	kmerlen_start -- integer, minimum length of k-mer in the list of svm weights
	kmerlen_end   -- integer, maximum length of k-mer in the list of sv weights

	Return:
	SVM score
	"""
	kmerlens = range(kmerlen_start, kmerlen_end+1)
	nkmerlens = len(kmerlens)

	score_norm_sum = 0

	for i in range(nkmerlens):
		score_norm = score_seq(s, svmwlist[i], kmerlens[i])
		score_norm_sum += score_norm
		
	return score_norm_sum


def main(argv = sys.argv):
	usage = "Usage: %prog [options] SVM_WEIGHTS TEST_SEQ"
	desc  = "1. take two files(one is in FASTA format to score, the other is SVM weight file generated from kmersvm_train.py) as input, 2. score each sequence in the given file"
	parser = optparse.OptionParser(usage=usage, description=desc)                                                                              
	parser.add_option("-o", dest="output", default="kmersvm_scores.out", \
  			help="set the name of output score file (default=kmersvm_scores.out)")

	parser.add_option("-q", dest="quiet", default=False, action="store_true", \
  			help="supress messages (default=false)")

	(options, args) = parser.parse_args()

	if len(args) == 0:
		parser.print_help()
		sys.exit(0)

	if len(args) != 2:
		parser.error("incorrect number of arguments")
		sys.exit(0)

	ktype_str = ["", "Spectrum", "Weighted Spectrums"]

	svmwf = args[0]
	seqf = args[1]

	seqs, sids = read_fastafile(seqf)
	svmwlist, params = read_svmwfile_wsk(svmwf)

	if options.quiet == False:
		sys.stderr.write('Options:\n')
		sys.stderr.write('  kernel-type: ' + str(params.kernel) + "." + ktype_str[params.kernel] + '\n')
		sys.stderr.write('  kmerlen: ' + str(params.kmerlen) + '\n')
		if params.kernel == 2:
			sys.stderr.write('  kmerlen2: ' + str(params.kmerlen2) + '\n')
		sys.stderr.write('  output: ' + options.output + '\n')
		sys.stderr.write('\n')

		sys.stderr.write('Input args:\n')
		sys.stderr.write('  SVM weights file: ' + svmwf + '\n')
		sys.stderr.write('  sequence file: ' + seqf + '\n')
		sys.stderr.write('\n')

		sys.stderr.write('numer of sequences to score: ' + str(len(seqs)) + '\n')
		sys.stderr.write('posteriorp A: ' + str(params.A) + '\n')
		sys.stderr.write('posteriorp B: ' + str(params.B) + '\n')
		sys.stderr.write('\n')

	f = open(options.output, 'w')
	f.write("\t".join(["#seq_id", "posterior_prob", "svm_score\n"]))

	kmerlen = params.kmerlen
	kmerlen2 = params.kmerlen2
	bias = params.bias
	A = params.A
	B = params.B
	for sidx in xrange(len(seqs)):
		s = seqs[sidx]
		score = score_seq_wsk(s, svmwlist, kmerlen, kmerlen2) + bias
		pp = 1/(1+numpy.exp(score*A+B))

		f.write("\t".join([ sids[sidx], str(pp), str(score)]) + "\n")

	f.close()

if __name__=='__main__': main()
