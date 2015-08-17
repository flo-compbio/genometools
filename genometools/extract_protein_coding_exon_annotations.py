#!/usr/bin/env python

# Copyright (c) 2015 Florian Wagner
#
# This file is part of GenomeTools.
#
# GenomeTools is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License, Version 3,
# as published by the Free Software Foundation.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

import sys
import os
import gzip
import csv
import re
import gzip
import argparse

def read_args_from_cmdline():
	parser = argparse.ArgumentParser(description='')

	parser.add_argument('-a','--annotation-file',default='-')
	parser.add_argument('-o','--output-file',required=True)
	parser.add_argument('-c','--chromosome-pattern',default=r'(?:\d\d?|MT|X|Y)$')
	parser.add_argument('-f','--field-name',default='exon')

	return parser.parse_args()

def open_plain_or_gzip(fn):
	try:
		gzip.open(fn).next()
		return gzip.open(fn)
	except IOError:
		return open(fn)

attr_sep = re.compile(r"(?<!\\)\s*;\s*") # use negative lookbehind to make sure we don't split on escaped semicolons ("\;")
def parse_attributes(s):
	''' parses the 9th field (attributes) of a GFF/GTF entry into a dictionary '''
	attr = {}
	atts = attr_sep.split(s)
	for a in atts:
		#print a
		kv = a.split(' ')
		if len(kv) == 2:
			k,v = kv
			v = v.strip('"')
			attr[k] = v
	return attr	

def main(args=None):

	if args is None:
		args = read_args_from_cmdline()

	input_file = args.annotation_file
	chrom_pat = re.compile(args.chromosome_pattern)
	field_name = args.field_name
	output_file = args.output_file

	print 'Regular expression used for filtering chromosome names:',chrom_pat.pattern

	chromosomes = set()
	excluded_chromosomes = set()
	i = 0
	exons = 0
	with open_plain_or_gzip(input_file) if input_file != '-' else sys.stdin as fh, open(output_file,'w') as ofh:
		#if i >= 500000: break
		reader = csv.reader(fh,dialect='excel-tab')
		writer = csv.writer(ofh,dialect='excel-tab',lineterminator=os.linesep,\
				quoting=csv.QUOTE_NONE,quotechar='|')
		for l in reader:
			i += 1
			if i % int(1e5) == 0:
				print '\r%d...' %(i), ; sys.stdout.flush() # report progress
			if len(l) > 1 and l[2] == field_name:
				attr = parse_attributes(l[8])
				type_ = attr['gene_biotype']
				if type_ in ['protein_coding','polymorphic_pseudogene']:

					# test whether chromosome is valid
					chrom = l[0]
					m = chrom_pat.match(chrom)
					if m is None:
						excluded_chromosomes.add(chrom)
						continue

					chromosomes.add(chrom)
					writer.writerow(l)
					exons += 1

	print "done (parsed %d lines)." %(i)

	print 
	print "Gene chromosomes (%d):" %(len(chromosomes))
	print "\t" + ', '.join(sorted(chromosomes))
	print
	print "Excluded chromosomes (%d):" %(len(excluded_chromosomes))
	print "\t" + ', '.join(sorted(excluded_chromosomes))
	print
	print "Total no. of exons:", exons

	return 0

if __name__ == '__main__':
	return_code = main()
	sys.exit(return_code)
