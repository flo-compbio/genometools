#!/usr/bin/env python

import sys
import argparse
import numpy as np
import gzip
import os

def read_args_from_cmdline():
	parser = argparse.ArgumentParser(description='Filter fasta.')
	parser.add_argument('-f','--fasta-file',default='-')
	parser.add_argument('-c','--chromosome-pattern',default=r'(?:\d\d?|MT|X|Y)$')
	parser.add_argumetn('-o','--output-file')
	args = parser.parse_args()
	return args

def open_plain_or_gzip(fn):
	try:
		gzip.open(fn).next()
		return gzip.open(fn)
	except IOError:
		return open(fn)

def main(args=None):

	if args is None:
		args = read_args_from_cmdline()

	input_file = args.fasta_file
	output_file = args.output_file
	chrom_pat = args.chromosome_pattern

	inside = False
	with open_plain_or_gzip(input_file) if input_file != '-' else sys.stdin as fh,\
			open(output_file,'w') if output_file != '-' else sys.stdout as ofh:

		for l in fh:
			if l[0] == '>':
				inside = True
				chrom = l[1:-1].split(' ',1)[0]
				m = chrom_pat.match(chrom)
				if m is None:
					inside = False # ignore this chromosome
			if inside:
				ofh.write(l)

	return 0

if __name__ == '__main__':
	return_code = main()
	sys.exit(return_code)
