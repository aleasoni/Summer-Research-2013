import os
import os.path
import sys
import random
import optparse

from bitstring import Bits, BitArray

def read_bed_file(filename):
	"""
	"""
	f = open(filename, 'r')
	lines = f.readlines()
	f.close()

	positions = []

	for line in lines:
		if line[0] == '#': 
			continue

		l = line.split()

		positions.append((l[0], int(l[1]), int(l[2])))

	return positions


def make_profile(positions, buildname, basedir):
	"""
	"""
	chrnames = sorted(set(map(lambda p: p[0], positions)))

	profiles = []
	for chrom in chrnames:
		idxf_gc = os.path.join(basedir, '.'.join([buildname, chrom, 'gc', 'out']))
		idxf_rpt = os.path.join(basedir, '.'.join([buildname, chrom, 'rpt', 'out']))

		#if os.path.exists(idxf_gc) == False or os.path.exists(idxf_rpt) == False:
		#	continue

		bits_gc = Bits(filename=idxf_gc)
		bits_rpt = Bits(filename=idxf_rpt)

		for pos in positions:
			if pos[0] != chrom:
				continue

			seqid = pos[0] + ':' + str(pos[1]+1) + '-' + str(pos[2])
			slen = pos[2]-pos[1]
			gc = bits_gc[pos[1]:pos[2]].count(True)
			rpt = bits_rpt[pos[1]:pos[2]].count(True)

			profiles.append((seqid, slen, gc, rpt))

	return profiles


def sample_sequences(positions, buildname, basedir, options):
	"""
	"""
	rpt_err = options.rpt_err
	gc_err = options.gc_err
	max_trys = options.max_trys
	norpt = options.norpt
	nogc = options.nogc

	chrnames = sorted(set(map(lambda p: p[0], positions)))
	profiles = make_profile(positions, buildname, basedir)

	excluded = []
	if options.skipfile:
		excluded = read_bed_file(options.skipfile, chr)

	f = open(options.output,"w")

	for chrom in chrnames:
		print chrom
		idxf_na = os.path.join(basedir, '.'.join([buildname, chrom, 'na', 'out']))
		idxf_gc = os.path.join(basedir, '.'.join([buildname, chrom, 'gc', 'out']))
		idxf_rpt = os.path.join(basedir, '.'.join([buildname, chrom, 'rpt', 'out']))

		bits_gc = Bits(filename=idxf_gc)
		bits_rpt = Bits(filename=idxf_rpt)

		#this bit array is used to mark positions that are excluded from sampling
		#this will be updated as we sample more sequences in order to prevent sampled sequences from overlapping
		bits_na = BitArray(filename=idxf_na)

		#mark excluded regions
		for pos in excluded:
			if pos[0] != chrom:
				continue
			bits_na.set(True, range(pos[1], pos[2]))
			npos+=1

		npos = 0
		#mark positive regions
		for pos in positions:
			if pos[0] != chrom:
				continue
			bits_na.set(True, range(pos[1], pos[2]))
			npos+=1

		if options.count == 0:
			count = options.fold*npos

		sampled_cnt = 0
		while sampled_cnt < count:
			sampled_prof = random.choice(profiles)
			sampled_len = sampled_prof[1]
			sampled_gc = sampled_prof[2]
			sampled_rpt = sampled_prof[3]

			rpt_err_allowed = int(rpt_err*sampled_len)
			gc_err_allowed = int(gc_err*sampled_len)
			trys = 0
			while trys < max_trys:
				trys += 1

				pos = random.randint(1, bits_na.length - sampled_len)
				pos_e = pos+sampled_len
		
				#if bits_na.any(True, range(pos, pos_e)):
				#	continue

				if bits_na[pos:pos_e].count(True) > 0:
					continue

				if not norpt:
					pos_rpt = bits_rpt[pos:pos_e].count(True)
					if abs(sampled_rpt - pos_rpt) > rpt_err_allowed:
						continue

				if not nogc:
					pos_gc = bits_gc[pos:pos_e].count(True)
					if abs(sampled_gc - pos_gc) > gc_err_allowed:
						continue

				#accept the sampled position

				#mark the sampled regions
				bits_na.set(True, range(pos, pos_e))

				f.write('\t'.join([chrom, str(pos), str(pos_e)]) + '\n')

				sampled_cnt += 1

				print trys, chrom, pos, pos_e, sampled_len, pos_rpt, sampled_rpt, pos_gc, sampled_gc
				break
			else:
				print "maximum trys reached"

	f.close()


def main(argv=sys.argv):
	usage = "usage: %prog [options] <Input Bed File> <Genome Build Name> <Base Directory>"

	desc  = "generate null sequences"
	parser = optparse.OptionParser(usage=usage, description=desc)

	parser.add_option("-x", dest="fold", type="int", \
		default = 1, help="number of sequence to sample, FOLD times of given dataset (default=1)")

	parser.add_option("-c", dest="count", type="int", \
		default=0, help="number of sequences to sample, override -x option (default=NA, obsolete)")

	parser.add_option("-r", dest="rseed", type="int", \
		default=1, help="random number seed (default=1)")

	parser.add_option("-g", dest="gc_err", type="float", \
		default=0.02, help="GC errors allowed (default=0.02)")

	parser.add_option("-t", dest="rpt_err", type="float", \
		default=0.02, help="RPT errors allowed (default=0.02)")

	parser.add_option("-s", dest="skipfile", \
		default="", help="filename that contains regions to be excluded (default=NA)")

	parser.add_option("-G", dest="nogc", action="store_true", \
		default=False, help="do not match gc-contents")

	parser.add_option("-R", dest="norpt", action="store_true", \
		default=False, help="do not match repeats")

	parser.add_option("-m", dest="max_trys", type="int", \
		default=10000, help="number of maximum trys to sample of one sequence (default=10000)")

	parser.add_option("-o", dest="output", default="null_seqs.bed", \
  			help="set the name of output file (default=null_seqs.bed)")

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


	posfile = args[0]
	buildname = args[1]
	basedir = args[2]

	random.seed(options.rseed)

	positions = read_bed_file(posfile)

	sample_sequences(positions, buildname, basedir, options)

		
if __name__ == "__main__": main()
