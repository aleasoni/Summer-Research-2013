#!/usr/bin/env python
"""
	kmersvm_train.py; train a support vector machine using shogun toolbox
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
import optparse
import random
import numpy
from math import log, exp

from libkmersvm import *
try:
	from shogun.PreProc import SortWordString, SortUlongString
except ImportError:
	from shogun.Preprocessor import SortWordString, SortUlongString
from shogun.Kernel import CommWordStringKernel, CommUlongStringKernel, \
		CombinedKernel
		
from shogun.Features import StringWordFeatures, StringUlongFeatures, \
		StringCharFeatures, CombinedFeatures, DNA, Labels
from shogun.Classifier import MSG_INFO, MSG_ERROR
try:
	from shogun.Classifier import SVMLight
except ImportError:
	from shogun.Classifier import LibSVM

"""
global variables
"""
g_kmers = []
g_rcmap = []


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


def non_redundant_word_features(feats, kmerlen):
	"""convert the features from Shogun toolbox to non-redundant word features (handle reverse complements)
	Arguments:
	feats -- StringWordFeatures
	kmerlen -- integer, length of k-mer

	Return:
	StringWordFeatures after converting reverse complement k-mer ids
	"""

	rcmap = g_rcmap

	for i in xrange(feats.get_num_vectors()):
		nf = [rcmap[int(kmerid)] for kmerid in feats.get_feature_vector(i)]

		feats.set_feature_vector(numpy.array(nf, numpy.dtype('u2')), i)

	preproc = SortWordString()
	preproc.init(feats)
	try:
		feats.add_preproc(preproc)
		feats.apply_preproc()
	except AttributeError:
		feats.add_preprocessor(preproc)
		feats.apply_preprocessor()	

	return feats


def non_redundant_ulong_features(feats, kmerlen):
	"""convert the features from Shogun toolbox to non-redundant ulong features
	Arguments:
	feats -- StringUlongFeatures
	kmerlen -- integer, length of k-mer

	Return:
	StringUlongFeatures after converting reverse complement k-mer ids
	"""

	for i in xrange(feats.get_num_vectors()):
		nf = [get_rcmap(int(kmerid), kmerlen) \
				for kmerid in feats.get_feature_vector(i)]

		feats.set_feature_vector(numpy.array(nf, numpy.dtype('u8')), i)

	preproc = SortUlongString()
	preproc.init(feats)
	try:
		feats.add_preproc(preproc)
		feats.apply_preproc()
	except AttributeError:
		feats.add_preprocessor(preproc)
		feats.apply_preprocessor()

	return feats


def svm_learn(kernel, labels, options):
	"""train SVM using SVMLight or LibSVM

	Arguments:
	kernel -- kernel object from Shogun toolbox
	lebels -- list of labels
	options -- object containing option data 

	Return:
	trained svm object 
	"""

	try: 
		svm=SVMLight(options.svmC, kernel, Labels(numpy.array(labels, dtype=numpy.double)))
	except NameError:
		svm=LibSVM(options.svmC, kernel, Labels(numpy.array(labels, dtype=numpy.double)))

	if options.quiet == False:
		svm.io.set_loglevel(MSG_INFO)
		svm.io.set_target_to_stderr()

	svm.set_epsilon(options.epsilon)
	svm.parallel.set_num_threads(1)
	if options.weight != 1.0:
		svm.set_C(options.svmC, options.svmC*options.weight)
	svm.train()

	if options.quiet == False:
		svm.io.set_loglevel(MSG_ERROR)

	return svm


def _get_spectrum_features(seqs, kmerlen):
	"""generate spectrum features (internal)

	Arguments:
	seqs -- list of sequences 
	kmerlen -- integer, length of k-mer

	Return:
	StringWord(Ulong)Features after treatment of redundant reverse complement k-mers
	"""

	char_feats = StringCharFeatures(seqs, DNA)

	if kmerlen <= 8:
		string_features = StringWordFeatures
		non_redundant_features = non_redundant_word_features
	else:
		string_features = StringUlongFeatures
		non_redundant_features = non_redundant_ulong_features
	
	feats = string_features(DNA)
	feats.obtain_from_char(char_feats, kmerlen-1, kmerlen, 0, False)
	return non_redundant_features(feats, kmerlen)


def get_spectrum_features(seqs, options):
	"""generate spectrum features (wrapper)
	"""
	return _get_spectrum_features(seqs, options.kmerlen)


def get_weighted_spectrum_features(seqs, options):
	"""generate weighted spectrum features
	"""
	global g_kmers
	global g_rcmap

	subfeats_list = []

	for k in xrange(options.kmerlen, options.kmerlen2+1):
		char_feats = StringCharFeatures(seqs, DNA)
		if k <= 8:
			g_kmers = generate_kmers(k)
			g_rcmap = generate_rcmap_table(k, g_kmers)

		subfeats = _get_spectrum_features(seqs, k)
		subfeats_list.append(subfeats)

	return subfeats_list


def get_spectrum_kernel(feats, options):
	"""build spectrum kernel with non-redundant k-mer list (removing reverse complement)

	Arguments:
	feats -- feature object
	options -- object containing option data 

	Return:
	StringWord(Ulong)Features, CommWord(Ulong)StringKernel
	"""
	if options.kmerlen <= 8:
		return CommWordStringKernel(feats, feats)
	else:
		return CommUlongStringKernel(feats, feats)


def get_weighted_spectrum_kernel(subfeats_list, options):
	"""build weighted spectrum kernel with non-redundant k-mer list (removing reverse complement)

	Arguments:
	subfeats_list -- list of sub-feature objects
	options -- object containing option data 

	Return:
	CombinedFeatures of StringWord(Ulong)Features, CombinedKernel of CommWord(Ulong)StringKernel 
	"""
	kmerlen = options.kmerlen
	kmerlen2 = options.kmerlen2

	subkernels = 0
	kernel = CombinedKernel()
	feats = CombinedFeatures()

	for subfeats in subfeats_list:
		feats.append_feature_obj(subfeats)

	for k in xrange(kmerlen, kmerlen2+1):
		if k <= 8:
			subkernel = CommWordStringKernel(10, False)
		else:
			subkernel = CommUlongStringKernel(10, False)

		kernel.append_kernel(subkernel)
		subkernels+=1

	kernel.init(feats, feats)

	kernel.set_subkernel_weights(numpy.array([1/float(subkernels)]*subkernels, numpy.dtype('float64')))

	return kernel
		

def init_spectrum_kernel(kern, feats_lhs, feats_rhs):
	"""initialize spectrum kernel (wrapper function)
	"""
	kern.init(feats_lhs, feats_rhs)


def init_weighted_spectrum_kernel(kern, subfeats_list_lhs, subfeats_list_rhs):
	"""initialize weighted spectrum kernel (wrapper function)
	"""
	feats_lhs = CombinedFeatures()
	feats_rhs = CombinedFeatures()

	for subfeats in subfeats_list_lhs:
		feats_lhs.append_feature_obj(subfeats)

	for subfeats in subfeats_list_rhs:
		feats_rhs.append_feature_obj(subfeats)

	kern.init(feats_lhs, feats_rhs)


def get_sksvm_weights(svm, feats, options):
	"""calculate the SVM weight vector of spectrum kernel
	"""
	kmerlen = options.kmerlen
	alphas = svm.get_alphas()
	support_vector_ids = svm.get_support_vectors()

	w = numpy.array([0]*(2**(2*kmerlen)), numpy.double)

	for i in xrange(len(alphas)):
		x = [0]*(2**(2*kmerlen))
		for kmerid in feats.get_feature_vector(int(support_vector_ids[i])):
			x[int(kmerid)] += 1
		x = numpy.array(x, numpy.double)
		w += (alphas[i]*x/numpy.sqrt(numpy.sum(x**2)))
	
	return w


def get_wsksvm_weights(svm, subfeats_list, options):
	"""calculate the SVM weight vector of weighted spectrum kernel
	"""
	kmerlen = options.kmerlen
	kmerlen2 = options.kmerlen2
	alphas = svm.get_alphas()
	support_vector_ids = svm.get_support_vectors()
	kmerlens = range(kmerlen, kmerlen2+1)

	weights = []
	for idx in xrange(len(kmerlens)):
		subfeats = subfeats_list[idx]

		k = kmerlens[idx]
		w = numpy.array([0]*(2**(2*k)), numpy.double)

		for i in xrange(len(alphas)):
			x = [0]*(2**(2*k))
			for kmerid in subfeats.get_feature_vector(int(support_vector_ids[i])):
				x[int(kmerid)] += 1
			x = numpy.array(x, numpy.double)
			w += (alphas[i]*x/numpy.sqrt(numpy.sum(x**2)))
	
		w /= len(kmerlens)
		weights.append(w)
	
	return weights 


def save_header(f, bias, A, B, options):
	f.write("#parameters:\n")
	f.write("#kernel=" + str(options.ktype) + "\n")
	f.write("#kmerlen=" + str(options.kmerlen) + "\n")
	if options.ktype == 2:
		f.write("#kmerlen2=" + str(options.kmerlen2) + "\n")
	f.write("#bias=" + str(bias) + "\n")
	f.write("#A=" + str(A) + "\n")
	f.write("#B=" + str(B) + "\n")
	f.write("#NOTE: k-mers with large negative weights are also important. They can be found at the bottom of the list.\n")
	f.write("#k-mer\trevcomp\tSVM-weight\n")


def save_sksvm_weights(w, bias, A, B, options):
	"""save the SVM weight vector from spectrum kernel
	"""
	output = options.outputname + "_weights.out"
	kmerlen = options.kmerlen

	f = open(output, 'w')
	save_header(f, bias, A, B, options)

	global g_kmers
	global g_rcmap

	if options.sort:
		w_sorted = sorted(zip(range(len(w)), w), key=lambda x: x[1], reverse=True)
	else:
		w_sorted = zip(range(len(w)), w)

	if kmerlen <= 8:
		for i in map(lambda x: x[0], w_sorted): 
			if i == g_rcmap[i]:
				f.write('\t'.join( [g_kmers[i], revcomp(g_kmers[i]), str(w[i])] ) + '\n')
	else:
		for i in map(lambda x: x[0], w_sorted): 
			if i == get_rcmap(i, kmerlen):
				kmer = kmerid2kmer(i, kmerlen)
				f.write('\t'.join( [kmer, revcomp(kmer), str(w[i])] ) + '\n')

	f.close()


def save_wsksvm_weights(w, bias, A, B, options):
	"""save the SVM weight vector from weighted spectrum kernel
	"""
	output = options.outputname + "_weights.out"
	kmerlen = options.kmerlen
	kmerlen2 = options.kmerlen2

	f = open(output, 'w')
	save_header(f, bias, A, B, options)

	global g_kmers
	global g_rcmap

	kmerlens = range(kmerlen, kmerlen2+1)
	for idx in xrange(len(kmerlens)):
		k = kmerlens[idx]
		subw = w[idx]

		if options.sort:
			subw_sorted = sorted(zip(range(len(subw)), subw), key=lambda x: x[1], reverse=True)
		else:
			subw_sorted = zip(range(len(subw)), subw)

		if k <= 8:
			g_kmers = generate_kmers(k)
			g_rcmap = generate_rcmap_table(k, g_kmers)
			for i in map(lambda x: x[0], subw_sorted): 
				if i == g_rcmap[i]:
					f.write('\t'.join( [g_kmers[i], revcomp(g_kmers[i]), str(subw[i])] ) + "\n")
		else:
			for i in map(lambda x: x[0], subw_sorted): 
				if i == get_rcmap(i, k):
					kmer = kmerid2kmer(i, k)
					f.write('\t'.join( [kmers, revcomp(kmers), str(subw[i])] ) + "\n")

	f.close()


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


def LMAI(svms, labels, prior0, prior1):
	"""fitting svms to sigmoid function (improved version introduced by Lin 2003)

	Arguments:
	svms -- list of svm scores
	labels -- list of labels
	prior0 -- prior of negative set
	prior1 -- prior of positive set

	Return:
	A, B parameter of 1/(1+exp(A*SVM+B))
	"""

	#parameter settings
	maxiter = 100
	minstep = 1e-10
	sigma = 1e-3

	hiTarget = (prior1+1.0)/float(prior1+2.0)
	loTarget = 1/float(prior0+2.0)

	t = [0]*len(labels)
	for i in xrange(len(labels)):
		if labels[i] == 1:
			t[i] = hiTarget
		else:
			t[i] = loTarget

	A = 0.0
	B = log((prior0+1.0)/float(prior1+1.0))
	fval = 0.0

	for i in xrange(len(labels)):
		fApB = svms[i]*A+B
		if fApB >= 0:
			fval += (t[i]*fApB+log(1+exp(-fApB)))
		else:
			fval += ((t[i]-1)*fApB+log(1+exp(fApB)))


	for it in xrange(maxiter):
		#print "iteration:", it
		#Update Graidient and Hessian (use H'= H + sigma I)
		h11 = sigma
		h22 = sigma
		h21 = 0.0
		g1 = 0.0
		g2 = 0.0

		for i in xrange(len(labels)):
			fApB = svms[i]*A+B
			if fApB >= 0:
				p = exp(-fApB) / float(1.0+exp(-fApB))
				q = 1.0 / float(1.0 + exp(-fApB))
			else:
				p = 1.0 / float(1.0 + exp(fApB))
				q = exp(fApB) / float(1.0+exp(fApB))
			d2 = p*q
			h11 += (svms[i]*svms[i]*d2)
			h22 += d2
			h21 += (svms[i]*d2)
			d1 = t[i]-p
			g1 += (svms[i]*d1)
			g2 += d1

		#Stopping criteria
		if (abs(g1)<1e-5) and (abs(g2)<1e-5):
			break

		det = h11*h22-h21*h21
		dA = -(h22*g1-h21*g2)/float(det)
		dB = -(-h21*g1+h11*g2)/float(det)
		gd = g1*dA+g2*dB
		stepsize=1
		while stepsize >= minstep:
			newA = A+stepsize*dA
			newB = B+stepsize*dB
			newf = 0.0

			for i in xrange(len(labels)):
				fApB = svms[i]*newA+newB
				if fApB >= 0:
					newf += (t[i]*fApB + log(1+exp(-fApB)))
				else:
					newf += ((t[i]-1)*fApB + log(1+exp(fApB)))

			if newf < (fval+0.0001*stepsize*gd):
				A=newA
				B=newB
				fval=newf
				break
			else:
				stepsize=stepsize/float(2.0)

		#Line search failes
		if stepsize < minstep:
			#print "Line search fails"
			break

	#if it >= maxiter:
	#	print "Reaching maximum iterations"

	return A, B


def wsksvm_classify(seqs, svm, kern, feats, options):
	feats_te = get_weighted_spectrum_features(seqs, options)
	init_weighted_spectrum_kernel(kern, feats, feats_te)

	return svm.apply().get_labels().tolist()


def score_seq(s, svmw, kmerlen):
	"""calculate SVM score of given sequence using single set of svm weights

	Arguments:
	s -- string, DNA sequence
	svmw -- numpy array, SVM weights 
	kmerlen -- integer, length of k-mer of SVM weight

	Return:
	SVM score
	"""

	global g_rcmap
	kmer2kmerid_func = kmer2kmerid

	x = [0]*(2**(2*kmerlen))
	for j in xrange(len(s)-kmerlen+1):
		x[ g_rcmap[kmer2kmerid_func(s[j:j+kmerlen], kmerlen)] ] += 1

	x = numpy.array(x, numpy.double)
	score_norm = numpy.dot(svmw, x)/numpy.sqrt(numpy.sum(x**2))

	return score_norm


def sksvm_classify(seqs, svm, kern, feats, options):
	"""classify the given sequences
	"""
	if options.kmerlen <= 8:
		#this is much faster when the length of kmer is short, and SVs are many
		svmw = get_sksvm_weights(svm, feats, options)
		return [score_seq(s, svmw, options.kmerlen)+svm.get_bias() for s in seqs]
	else:
		feats_te = get_spectrum_features(seqs, options)
		init_spectrum_kernel(kern, feats, feats_te)

		return svm.apply().get_labels().tolist()


def main(argv = sys.argv):
	usage = "Usage: %prog [options] POSITIVE_SEQ NEGATIVE_SEQ"
	desc  = "1. take two files(FASTA format) as input, 2. train an SVM and store the trained SVM weights"
	parser = optparse.OptionParser(usage=usage, description=desc)
	parser.add_option("-t", dest="ktype", type="int", default=1, \
			help="set the type of kernel, 1:Spectrum, 2:Weighted Spectrums (default=1.Spectrum)")

	parser.add_option("-C", dest="svmC", type="float", default=1, \
			help="set the regularization parameter svmC (default=1)")

	parser.add_option("-e", dest="epsilon", type="float", default=0.00001, \
			help="set the precision parameter epsilon (default=0.00001)")

	parser.add_option("-w", dest="weight", type="float", default=0.0, \
			help="set the weight for positive set (default=auto, 1+log(N/P))")

	parser.add_option("-k", dest="kmerlen", type="int",default=6, \
			help="set the (min) length of k-mer for (weighted) spectrum kernel (default = 6)")

	parser.add_option("-K", dest="kmerlen2", type="int",default=8, \
			help="set the max length of k-mer for weighted spectrum kernel (default = 8)")

	parser.add_option("-n", dest="outputname", default="kmersvm_output", \
  			help="set the name of output files (default=kmersvm_output)")

	parser.add_option("-v", dest="ncv", type="int", default=0, \
			help="if set, it will perform N-fold cross-validation and generate a prediction file (default = 0)")

	parser.add_option("-p", dest="posteriorp", default=False, action="store_true", \
  			help="estimate parameters for posterior probability with N-CV. this option requires -v option to be set (default=false)")

	parser.add_option("-r", dest="rseed", type="int", default=1, \
			help="set the random number seed for cross-validation (-p option) (default=1)")

	parser.add_option("-q", dest="quiet", default=False, action="store_true", \
  			help="supress messages (default=false)")

	parser.add_option("-s", dest="sort", default=False, action="store_true", \
  			help="sort the kmers by absolute values of SVM weights (default=false)")

	ktype_str = ["", "Spectrum", "Weighted Spectrums"]

	(options, args) = parser.parse_args()

	if len(args) == 0:
		parser.print_help()
		sys.exit(0)

	if len(args) != 2:
		parser.error("incorrect number of arguments")
		parser.print_help()
		sys.exit(0)

	if options.posteriorp and options.ncv == 0:
		parser.error("posterior probability estimation requires N-fold CV process (-v option should be set)")
		parser.print_help()
		sys.exit(0)

	random.seed(options.rseed)

	"""
	set global variable
	"""
	if (options.ktype == 1) and (options.kmerlen <= 8):
		global g_kmers
		global g_rcmap

		g_kmers = generate_kmers(options.kmerlen)
		g_rcmap = generate_rcmap_table(options.kmerlen, g_kmers)
	
	posf = args[0]
	negf = args[1]
	
	seqs_pos, sids_pos = read_fastafile(posf)
	seqs_neg, sids_neg = read_fastafile(negf)
	npos = len(seqs_pos)
	nneg = len(seqs_neg)
	seqs = seqs_pos + seqs_neg
	sids = sids_pos + sids_neg

	if options.weight == 0:
		options.weight = 1 + log(nneg/npos)

	if options.quiet == False:
		sys.stderr.write('SVM parameters:\n')
		sys.stderr.write('  kernel-type: ' + str(options.ktype) + "." + ktype_str[options.ktype] + '\n')
		sys.stderr.write('  svm-C: ' + str(options.svmC) + '\n')
		sys.stderr.write('  epsilon: ' + str(options.epsilon) + '\n')
		sys.stderr.write('  weight: ' + str(options.weight) + '\n')
		sys.stderr.write('\n')

		sys.stderr.write('Other options:\n')
		sys.stderr.write('  kmerlen: ' + str(options.kmerlen) + '\n')
		if options.ktype == 2:
			sys.stderr.write('  kmerlen2: ' + str(options.kmerlen2) + '\n')
		sys.stderr.write('  outputname: ' + options.outputname + '\n')
		sys.stderr.write('  posteriorp: ' + str(options.posteriorp) + '\n')
		if options.ncv > 0:
			sys.stderr.write('  ncv: ' + str(options.ncv) + '\n')
			sys.stderr.write('  rseed: ' + str(options.rseed) + '\n')
		sys.stderr.write('  sorted-weight: ' + str(options.sort) + '\n')

		sys.stderr.write('\n')

		sys.stderr.write('Input args:\n')
		sys.stderr.write('  positive sequence file: ' + posf + '\n')
		sys.stderr.write('  negative sequence file: ' + negf + '\n')
		sys.stderr.write('\n')

		sys.stderr.write('numer of total positive seqs: ' + str(npos) + '\n')
		sys.stderr.write('numer of total negative seqs: ' + str(nneg) + '\n')
		sys.stderr.write('\n')

	#generate labels
	labels = [1]*npos + [-1]*nneg

	if options.ktype == 1:
		get_features = get_spectrum_features
		get_kernel = get_spectrum_kernel
		get_weights = get_sksvm_weights
		save_weights = save_sksvm_weights
		svm_classify = sksvm_classify
	elif options.ktype == 2:
		get_features = get_weighted_spectrum_features
		get_kernel = get_weighted_spectrum_kernel
		get_weights = get_wsksvm_weights
		save_weights = save_wsksvm_weights
		svm_classify = wsksvm_classify
	else:
		sys.stderr.write('..unknown kernel..\n')
		sys.exit(0)

	A = B = 0
	if options.ncv > 0:
		if options.quiet == False:
			sys.stderr.write('..Cross-validation\n')

		cvlist = generate_cv_list(options.ncv, npos, nneg)
		labels_cv = []
		preds_cv = []
		sids_cv = []
		indices_cv = []
		for icv in xrange(options.ncv):
			#split data into training and test set
			seqs_tr, seqs_te = split_cv_list(cvlist, icv, seqs) 
			labs_tr, labs_te = split_cv_list(cvlist, icv, labels)
			sids_tr, sids_te = split_cv_list(cvlist, icv, sids)
			indices_tr, indices_te = split_cv_list(cvlist, icv, range(len(seqs)))

			#train SVM
			feats_tr = get_features(seqs_tr, options)
			kernel_tr = get_kernel(feats_tr, options)
			svm_cv = svm_learn(kernel_tr, labs_tr, options)

			preds_cv = preds_cv + svm_classify(seqs_te, svm_cv, kernel_tr, feats_tr, options)
			
			labels_cv = labels_cv + labs_te
			sids_cv = sids_cv + sids_te
			indices_cv = indices_cv + indices_te

		output_cvpred = options.outputname + "_cvpred.out"
		prediction_results = sorted(zip(indices_cv, sids_cv, preds_cv, labels_cv), key=lambda p: p[0])
		save_predictions(output_cvpred, prediction_results, cvlist)

		if options.posteriorp:
			A, B = LMAI(preds_cv, labels_cv, labels_cv.count(-1), labels_cv.count(1))

			if options.quiet == False:
				sys.stderr.write('Estimated Parameters:\n')
				sys.stderr.write('  A: ' + str(A) + '\n')
				sys.stderr.write('  B: ' + str(B) + '\n')

	if options.quiet == False:
		sys.stderr.write('..SVM weights\n')

	feats = get_features(seqs, options)
	kernel = get_kernel(feats, options)
	svm = svm_learn(kernel, labels, options)
	w = get_weights(svm, feats, options)
	b = svm.get_bias()

	save_weights(w, b, A, B, options)

if __name__=='__main__': main()
