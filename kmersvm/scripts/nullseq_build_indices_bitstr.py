import os
import sys
import tarfile

from bitstring import Bits


def clear_indexes(sid, buildname):
	na = '.'.join([buildname, sid, "na", "out"])
	gc = '.'.join([buildname, sid, "gc", "out"])
	rpt = '.'.join([buildname, sid, "rpt", "out"])

	#truncate files
	for fn in (na, gc, rpt):
		f = open(fn, 'wb')
		f.close()


def append_indexes(seq, sid, buildname):
	na = '.'.join([buildname, sid, "na", "out"])
	gc = '.'.join([buildname, sid, "gc", "out"])
	rpt = '.'.join([buildname, sid, "rpt", "out"])

	f = open(na, 'ab')
	Bits(map(lambda c: c in 'N', seq)).tofile(f)
	f.close()

	f = open(gc, 'ab')
	Bits(map(lambda c: c in 'cgCG', seq)).tofile(f)
	f.close()

	f = open(rpt, 'ab')
	Bits(map(lambda c: c in 'acgt', seq)).tofile(f)
	f.close()


def build_indexes(filename, buildname):
	save_interval = 8*16*1024

	try:
		f = open(filename, 'r')

		seq = [] 
		sid = ''
		nlines = 0
		for line in f:
			if line[0] == '>':
				if sid: 
					append_indexes("".join(seq), sid, buildname)
					seq = []

				sid = line[1:].rstrip('\n').split()[0]
				clear_indexes(sid, buildname)
			else:
				nlines += 1
				seq.append(line.rstrip('\n'))

				if nlines % save_interval == 0:
					append_indexes("".join(seq), sid, buildname)
					seq = []

		#the last remaining sequence 
		append_indexes("".join(seq), sid, buildname)

	except IOError, (errno, strerror):
		print "I/O error(%d): %s" % (errno, strerror)
		sys.exit(0)


def main(argv=sys.argv):

	chrom_tar_file = argv[1]
	genome = argv[2]

	tar = tarfile.open(chrom_tar_file)

	for f in tar.getnames():
		print "processing", f
		tar.extract(f)
		build_indexes(f, genome)
		os.remove(f)

if __name__ == "__main__": main()
