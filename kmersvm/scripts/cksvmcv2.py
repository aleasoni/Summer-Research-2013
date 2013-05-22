import sys
import random
import numpy
import optparse
import itertools
from math import log, exp

from shogun.Kernel import CustomKernel, SqrtDiagKernelNormalizer
from shogun.Features import Labels
from shogun.Classifier import SVMLight, LibSVM, MSG_INFO, MSG_ERROR, MSG_DEBUG

#global variable
g_svm_bias = 0

def read_fastafile(filename, subs=True):
	"""Read a file in FASTA format

	Arguments:
	filename -- string, the name of the sequence file in FASTA format

	Return: 
	list of sequences, list of sequence ids
	"""
	id = '' 
	seqids = []
	seqs = []

	try:
		f = open(filename, 'r')
		lines = f.readlines()
		f.close()

	except IOError, (errno, strerror):
		print "I/O error(%d): %s" % (errno, strerror)
		sys.exit(0)

	seq = [] 
	for line in lines:
		if line[0] == '>':
			seqid = '_'.join(line[1:].rstrip('\n').split())
			seqids.append(seqid)
			if seq != []: seqs.append("".join(seq))
			seq = []
		else:
			if subs:
				seq.append(line.rstrip('\n').replace('N', 'A').upper())
			else:
				seq.append(line.rstrip('\n').upper())

	if seq != []:
		seqs.append("".join(seq))

	return seqs, seqids


def split_cv_list(cvlist, icv, data):
	"""
	"""
	tr_data = []
	te_data = []

	for i in range(len(data)):
		if cvlist[i] == icv:
			te_data.append(data[i])
		else:
			tr_data.append(data[i])
	
	return tr_data, te_data


def generate_cv_list(ncv, n1, n2):
	"""
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


def get_cv_list(sids, filename):
	try:
		f = open(filename, 'r')
		lines = f.readlines()
		f.close()

	except IOError, (errno, strerror):
		print "I/O error(%d): %s" % (errno, strerror)
		sys.exit(0)

	cv = []
	cvsids = []
	if len(sids) > 0:
		id2icv = {}
		for line in lines:
			f = line.rstrip('\n').split()
			id2icv[f[0]] = int(f[1])

		cv = [0] * len(sids)

		for idx in range(len(sids)):
			cvsids.append(sids[idx])
			cv[idx] = id2icv[sids[idx]]
	else:
		#assume that cvfile has the same order as matrix if sid is not provided
		cv = []
		for line in lines:
			f = line.rstrip('\n').split()
			cvsids.append(f[0])
			cv.append(int(f[1]))

	return cv, cvsids


def get_training_triangle_matrix (mat, cvlist, icv):
	idx_tr = []

	for i in xrange(len(mat)):
		if cvlist[i] != icv:
			idx_tr.append(i)

	#tr_triu_indices1 = [x for x in idx_tr for y in idx_tr if y>=x]
	#tr_triu_indices2 = [y for x in idx_tr for y in idx_tr if y>=x]

	ntr = len(idx_tr)
	idx_tr = numpy.array(idx_tr, dtype=numpy.int)
	rpttimes = list(range(1, ntr+1))
	rpttimes.reverse()
	tr_triu_indices1 = numpy.repeat(idx_tr, rpttimes)
	tr_triu_indices2 = numpy.empty(ntr*(ntr+1)/2, dtype=numpy.int)
	idx = 0
	for i in reversed(xrange(1, ntr+1)):
		tr_triu_indices2[idx:(idx+i)] = idx_tr[ntr-i:]
		idx += i

	#for i in xrange(1000000):
	#	print tr_triu_indices1a[i], tr_triu_indices1[i], tr_triu_indices2a[i], tr_triu_indices2[i]
	#sys.exit(0)

	#it should be upper triangle (the example was wrong!)
	return mat[(tr_triu_indices1, tr_triu_indices2)]


def get_training_full_matrix (mat, cvlist, icv):
	idx_tr = []

	for i in xrange(len(mat)):
		if cvlist[i] != icv:
			idx_tr.append(i)

	return mat[idx_tr, :][:, idx_tr]
	#return numpy.array([mat[(x,y)] for x in idx_tr for y in idx_tr]).reshape((len(idx_tr), len(idx_tr)))


def get_test_full_matrix (mat, cvlist, icv):
	colidx = []
	rowidx = []

	for i in xrange(len(mat)):
		if cvlist[i] != icv:
			rowidx.append(i)

	for i in xrange(len(mat)):
		if cvlist[i] == icv:
			colidx.append(i)

	
	return mat[rowidx, :][:, colidx]
	#return numpy.array([mat[(x,y)] for x in rowidx for y in colidx]).reshape((len(rowidx), len(colidx)))


def get_full_matrix(matrix_file, npos, nneg):
	#MODIFIED by dlee - 6.14.12
	#Mahmoud changed the output format such that it only contains lower triangle and empty upper (no 0s)
	#
	#input is lower triangle, but we need upper triangle HAHAHA, oh well, i am just going to make symmat
	ltmat = numpy.fromfile(matrix_file, dtype=numpy.float32, sep="\t")

	#MODIFIED by dlee - 9.22.12 
	#to speed up
	#ltmat_idx = 0
	#for row in xrange(npos+nneg):
	#	for col in xrange(row+1):
	#		mat[row, col] = ltmat[ltmat_idx]
	#		ltmat_idx += 1
	#ltmat = []

	mat = numpy.zeros((npos+nneg, npos+nneg))
	#tril_indices = ((x,y) for x in xrange(npos+nneg) for y in xrange(x+1))
	#for (row, col), val in itertools.izip(tril_indices, ltmat):
	#	mat[(row, col)] = val

	#unfortunately, memory error occurred in python 2.4
	#tril_indices1 = [x for x in indices for y in indices if y<=x]
	#tril_indices2 = [y for x in indices for y in indices if y<=x]
	tril_indices1 = numpy.repeat(xrange(npos+nneg), xrange(1,npos+nneg+1))
	tril_indices2 = numpy.empty((npos+nneg)*(npos+nneg+1)/2, dtype=numpy.int)
	idx = 0
	for i in xrange(1, npos+nneg+1):
		tril_indices2[idx:(idx+i)] = numpy.arange(i, dtype=numpy.int)
		idx += i

	sys.stderr.write(str(len(tril_indices1)) + '\n')
	sys.stderr.write(str(len(tril_indices2)) + '\n')
	sys.stderr.write(str(len(ltmat)) + '\n')
	mat[(tril_indices1, tril_indices2)] = ltmat

	#make symmetric mat
	mat.T[(tril_indices1, tril_indices2)] = ltmat

	return mat


def svm_learn(kernel, labels, svmC, epsilon, weight):
	"""
	"""
	try: 
		svm=SVMLight(svmC, kernel, Labels(numpy.array(labels, dtype=numpy.double)))
	except NameError:
		print 'No support for SVMLight available.'
		return

	svm.io.set_loglevel(MSG_INFO)
	svm.io.set_target_to_stderr()

	svm.set_epsilon(epsilon)
	svm.parallel.set_num_threads(1)
	if weight != 1.0:
		svm.set_C(svmC, svmC*weight)
	svm.train()
	svm.io.set_loglevel(MSG_ERROR)

	return svm


def save_cv_taining_result(svm, options, seqids, icv):
	sys.stderr.write('..writing outputs..\n')
	alphas = svm.get_alphas()
	support_vector_ids = svm.get_support_vectors()

	alphaf = '.'.join([options.alphaprefix, str(icv), "out"])

	f = open(alphaf, 'w')
	for i in xrange(len(alphas)):
		f.write(seqids[int(support_vector_ids[i])] + "\t")
		f.write(str(alphas[i]) + "\n")
	f.close()


def svm_train_classify(trimat_tr, seqids_tr, labels_tr, fullmat_te, options, icv):

	sys.stderr.write('..kernel building..\n')

	kernel=CustomKernel()
	kernel.set_triangle_kernel_matrix_from_triangle(trimat_tr)

	sys.stderr.write('..svm learning..\n')
	svm = svm_learn(kernel, labels_tr, options.svmC, options.epsilon, options.weight)

	global g_svm_bias
	g_svm_bias = svm.get_bias()

	if (options.alphaprefix != "") and (len(seqids_tr) > 0):
		save_cv_taining_result(svm, options, seqids_tr, icv)

	sys.stderr.write('..svm classifying..\n')
	kernel.set_full_kernel_matrix_from_full(fullmat_te)

	###################################################
	#for testing
	#alphas = svm.get_alphas()
	#svids  = svm.get_support_vectors()

	#for j in xrange(len(preds)):
	#	p = svm.get_bias()
	#	for i in xrange(len(alphas)):
	#		p += (alphas[i]*fullmat_te[int(svids[i]),j])

	#	print preds[j], p

	#sys.exit(0)
	###################################################

	return svm.classify().get_labels().tolist()


def LMAI(svms, labels, prior0, prior1):
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
			sys.stderr.write("Line search fails\n")
			break

	if it >= maxiter:
		sys.stderr.write("Reaching maximum iterations\n")

	return A, B


def write_pred_results(fh, ids, preds, labels, icv):
	#positives
	for i in range(len(preds)):
		if labels[i] == 1:
			fh.write("%s\t%.6f\t%d\t%d\n" % (ids[i], preds[i], labels[i], icv))

	#negatives
	for i in range(len(preds)):
		if labels[i] == -1:
			fh.write("%s\t%.6f\t%d\t%d\n" % (ids[i], preds[i], labels[i], icv))

def get_auc(preds, labs):
	preds_sig = map(lambda x: round(x, 4), preds)

	pred_tuples = zip(preds_sig, labs)
	pred_tuples_sorted = sorted(pred_tuples, key=lambda p: p[0], reverse=True)

	P = labs.count(1)
	N = labs.count(-1)

	auc = 0
	TP_prev = 0
	TP_curr = 0
	FP_prev = 0
	FP_curr = 0
	cutoff = pred_tuples_sorted[0][0]

	for i in xrange(len(pred_tuples_sorted)):
		if cutoff != pred_tuples_sorted[i][0]:
			auc += (TP_curr+TP_prev)*(FP_curr-FP_prev)
			TP_prev = TP_curr
			FP_prev = FP_curr
			cutoff = pred_tuples_sorted[i][0]

		if pred_tuples_sorted[i][1] == 1:
			TP_curr += 1
		else:
			FP_curr += 1
	
	if (TP_prev != TP_curr) or (FP_prev != FP_curr):
		auc += (TP_curr+TP_prev)*(FP_curr-FP_prev)

	return auc/float((2*P*N))

def main(argv=sys.argv) :
	usage = "Usage: %prog [options] KERNEL_MATRIX POSSEQF NEGSEQF OUTPUTF"
	desc  = "1. take a kernel matrix and the corresponding sequence files as input, 2. train SVM classifiers with n-fold cross-validation 3. classify and print the SVM prediction scores for the held-out test set to OUTPUT"

	parser = optparse.OptionParser(usage=usage, description=desc)                                                                              
	parser.add_option("-C", dest="svmC", type="float", default=1, \
			help="set the regularization parameter svmC (default=1)")

	parser.add_option("-e", dest="epsilon", type="float", default=0.00001, \
			help="set the precision parameter epsilon (default=0.00001)")

	parser.add_option("-w", dest="weight", type="float", default=1.0, \
			help="set the weight for positive set (default=1.0)")

	parser.add_option("-v", dest="ncv", type="int", default=5, \
			help="set the number of of cross-validation (default = 5)")

	parser.add_option("-i", dest="icv", type="int", \
			help="set a specific cross-validation set to run (from 0 to (NCV-1)). If not set, it will iterate all sets")

	parser.add_option("-r", dest="rseed", type="int", default=1, \
			help="set the random number seed (default=1)")

	parser.add_option("-p", dest="posteriorp", default=False, action="store_true",\
			help="if set, svm scores will be converted to posterior probabilities (default=False)")

	parser.add_option("-P", dest="paramprefix", default="",\
			help="set the prefix of output of the sigmoid parameter of each cv-training. this should be used with -p option. output will be ${paramprefix}.${icv}.out(default=NA)")

	parser.add_option("-a", dest="alphaprefix", default="",\
			help="set the prefix of output of the alpha file of cv-training (default=NA)")

	parser.add_option("-f", dest="cvfile", default="",\
			help="set the cv-file (if set, rseed will not be used)")

	(options, args) = parser.parse_args()

	if len(args) == 0:
		parser.print_help()
		sys.exit(0)

	if len(args) != 4:
		parser.error("incorrect number of arguments")
		sys.exit(0)

	matrix_file = args[0]
	posseqf = args[1]
	negseqf = args[2]
	outputf = args[3]

	sys.stderr.write('loading sequence files..\n')
	posseqs, posseqids = read_fastafile(posseqf)
	negseqs, negseqids = read_fastafile(negseqf)

	npos = len(posseqids)
	nneg = len(negseqids)

	seqids = posseqids + negseqids
	seqs = posseqs + negseqs
	labels = [1]*npos + [-1]*nneg

	random.seed(options.rseed)

	cvlist = []

	sys.stderr.write('loading kernel matrix..\n')

	#MOTIFIED 9.22.12
	mat = get_full_matrix(matrix_file, npos, nneg)

	#generate cross-validation list
	if options.cvfile == "":
		sys.stderr.write('generating cvlist..\n')
		cvlist = generate_cv_list(options.ncv, npos, nneg)
	else:
		sys.stderr.write('loading cvfile.. ' + options.cvfile + '\n')
		cvlist, cvsids = get_cv_list(seqids, options.cvfile)
		if len(seqids) == 0:
			seqids = cvsids

	fh = open(outputf, 'w')

	for icv in xrange(options.ncv):
		if (options.icv != None) and (icv != options.icv):
			continue

		sys.stderr.write('cross-validation: ' + str(icv) + '\n')
		labels_tr, labels_te = split_cv_list(cvlist, icv, labels)
		if len(seqids) > 0:
			seqids_tr, seqids_te = split_cv_list(cvlist, icv, seqids)
		else:
			seqids_tr = seqids_te = []

		#obtaining training matrix
		sys.stderr.write('..get training matrix..\n')
		trimat_tr = get_training_triangle_matrix(mat, cvlist, icv)
		fullmat_te = get_test_full_matrix(mat, cvlist, icv)

		preds_te = svm_train_classify(trimat_tr, seqids_tr, labels_tr, fullmat_te, options, icv)

		preds_idx = []
		for i in xrange(len(mat)):
			if cvlist[i] == icv:
				preds_idx.append(i)
		pps = preds_te

		if options.posteriorp == True:
			sys.stderr.write('...estimating sigmoid parameters using nested cross-validations...\n')

			cvlist2 = generate_cv_list(options.ncv, labels_tr.count(1), labels_tr.count(-1))
			fullmat_tr =  get_training_full_matrix(mat, cvlist, icv)

			preds_tr_te_all = []
			labels_tr_te_all = []
			for icv2 in range(options.ncv):
				labels_tr_tr, labels_tr_te = split_cv_list(cvlist2, icv2, labels_tr)
				trimat_tr_tr = get_training_triangle_matrix(fullmat_tr, cvlist2, icv2)
				fullmat_tr_te = get_test_full_matrix(fullmat_tr, cvlist2, icv2)

				preds_tr_te_all += svm_train_classify(trimat_tr_tr, [], labels_tr_tr, fullmat_tr_te, options, icv2)
				labels_tr_te_all += labels_tr_te

			AI, BI = LMAI(preds_tr_te_all, labels_tr_te_all, labels_tr_te_all.count(-1), labels_tr_te_all.count(1))

			if options.paramprefix != "":
				global g_svm_bias
				f = open('.'.join([options.paramprefix, str(icv), "out"]), 'w')
				f.write('A:\t' + str(AI) + '\n')
				f.write('B:\t' + str(BI) + '\n')
				f.write('bias:\t' + str(g_svm_bias) + '\n')
				f.close()

			sys.stderr.write('... A: ' + str(AI) + ', B: ' + str(BI) + "\n")
			pps = 1/(1+numpy.exp(numpy.array(preds_te, numpy.dtype('float64'))*AI+BI))

		if len(seqids) == 0:
			write_pred_results(fh, preds_idx, pps, labels_te, icv)
		else:
			write_pred_results(fh, seqids_te, pps, labels_te, icv)

		auc = get_auc(preds_te, labels_te)
		sys.stderr.write('auc = ' + str(auc) + '\n')

	fh.close()

if __name__ == '__main__': main()

