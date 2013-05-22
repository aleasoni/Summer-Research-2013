import os
import os.path
import sys
import optparse
import math
import re
from libkmersvm import *
		
def split(bed_file,options):
	split_f = open(options.output, 'w')	
	incr = options.incr
	size = options.size
	file = open(bed_file, 'rb')
	
	for line in file:		
		(name,start,length) = line.split('\t')
		start = int(start)
		length = int(length)
		end = size  + start

		while True:				
			coords = "".join([name,"\t",str(start),"\t",str(end),"\n"])
			split_f.write(coords)
			if end + incr >= length:
				end += incr-((end+incr)-length)
				start += incr	
				coords = "".join([name,"\t",str(start),"\t",str(end),"\n"])
				split_f.write(coords)	
				break						
			else:
				start += incr
				end += incr

				
def main(argv=sys.argv):
	usage = "usage: %prog <bed_file>"
	parser = optparse.OptionParser(usage=usage)

	parser.add_option("-s", dest="size", type="int", \
		default=1000, help="set chunk size")
	parser.add_option("-i", dest="incr", type="int", \
		default=500, help="set overlap size")
	parser.add_option("-o", dest="output", default="split_genome_output.bed", \
		help="output BED file (default is split_genome_output.bed)")
		
	(options, args) = parser.parse_args()
        if len(args) == 0:
		parser.print_help()
		sys.exit(0)
	
	bed_file = args[0]
	
	split(bed_file, options)	

if __name__ == "__main__": main()
