#!/usr/bin/env python

import sys
import os
import argparse
import gzip
import csv
import re

import pysam

def read_args_from_cmdline():
	parser = argparse.ArgumentParser(description='')
	parser.add_argument('-g','--genome-file',required=True)
	parser.add_argument('-c','--chromosome-pattern',default=r'(?:\d\d?|MT|X|Y)$')
	parser.add_argument('-o','--output-file',required=True)
	return parser.parse_args()

def main(args=None):

	if args is None:
		args = read_args_from_cmdline()

	genome_file = args.genome_file
	output_file = args.output_file
	chrom_pat = re.compile(args.chromosome_pattern)

	genome = pysam.Fastafile(genome_file)
	#chromosomes = set(args.chromosomes)
	#if chromosomes:
	#	print "Selected chromosomes:", ' '.join(sorted(chromosomes))

	ignored_chromosomes = set()
	chromlen = {}
	for chrom,length in zip(genome.references,genome.lengths):
		#print r,l
		m = chrom_pat.match(chrom)
		if m is None:
			ignored_chromosomes.add(chrom)
			continue

		chromlen[chrom] = length

	if ignored_chromosomes:
		print "Ignored chromosomes (%d): " %(len(ignored_chromosomes))
		print ', '.join(sorted(ignored_chromosomes))

	n = len(chromlen)
	print "Extracted chromosome lengths for %d chromosomes." %(n)

	with open(output_file,'w') as ofh:
		writer = csv.writer(ofh,dialect='excel-tab',lineterminator=os.linesep,quoting=csv.QUOTE_NONE)
		for chrom in sorted(chromlen.keys()):
			writer.writerow([chrom,str(chromlen[chrom])])

	return 0

if __name__ == '__main__':
	return_code = main()
	sys.exit(return_code)
