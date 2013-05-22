#!/usr/bin/env python2.7

import os
import os.path
import sys
import random
import optparse

from bitarray import bitarray

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


def bitarray_fromfile(filename):
	"""
	"""
	fh = open(filename, 'rb')
	bits = bitarray()
	bits.fromfile(fh)

	return bits, fh

def get_seqid(buildname, pos):
	return '_'.join( [buildname, pos[0], str(pos[1]), str(pos[2]), '+'] )

def make_profile(positions, buildname, basedir):
	"""
	"""
	chrnames = sorted(set(map(lambda p: p[0], positions)))

	profiles = {}
	for chrom in chrnames:
		idxf_gc = os.path.join(basedir, '.'.join([buildname, chrom, 'gc', 'out']))
		idxf_rpt = os.path.join(basedir, '.'.join([buildname, chrom, 'rpt', 'out']))

		#if os.path.exists(idxf_gc) == False or os.path.exists(idxf_rpt) == False:
		#	continue
		bits_gc, gcf = bitarray_fromfile(idxf_gc)
		bits_rpt, rptf = bitarray_fromfile(idxf_rpt)

		for pos in positions:
			if pos[0] != chrom:
				continue

			seqid = get_seqid(buildname, pos)
			slen = pos[2]-pos[1]
			gc = bits_gc[pos[1]:pos[2]].count(True)
			rpt = bits_rpt[pos[1]:pos[2]].count(True)

			profiles[seqid] = (slen, gc, rpt)

		gcf.close()
		rptf.close()

	return profiles

def main():
	parser = optparse.OptionParser()
	if len(sys.argv) != 5:
		print "Usage:", sys.argv[0], "BEDFILE BUILDNAME BASE_DIR OUT_FILE"
		sys.exit()
	parser.add_option("-o", dest="output", default="seq_profile.txt")
	(options,args) = parser.parse_args()
	
	bedfile = sys.argv[1]
	buildname = sys.argv[2]
	basedir = sys.argv[3]
	output = options.output

	positions = read_bed_file(bedfile)
	seqids = []
	for pos in positions:
		seqids.append(get_seqid(buildname, pos))

	profiles = make_profile(positions, buildname, basedir)

	f = open(output, 'w')
	f.write("\t".join(["#seq_id", "length", "GC content", "repeat fraction"]) + '\n')
	for seqid in seqids:
		prof = profiles[seqid]
		f.write('\t'.join( map(str, [seqid, prof[0], prof[1]/float(prof[0]), prof[2]/float(prof[0])]) ) + '\n')
	f.close()

if __name__ == "__main__": main()
