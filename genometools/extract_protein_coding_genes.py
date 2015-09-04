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

from collections import Counter

def read_args_from_cmdline():
	parser = argparse.ArgumentParser(description='')

	parser.add_argument('-a','--annotation-file',default='-')
	parser.add_argument('-o','--output-file',required=True)
	parser.add_argument('-c','--chromosome-pattern',default=r'(?:\d\d?|MT|X|Y)$')
	parser.add_argument('-f','--field-name',default='gene')

	#parser.add_argument('-e','--exclude-chromosomes',default=[],nargs='+')

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

	#if exclude_chromosomes:
	#	print "Excluding chromosomes %s..." %(', '.join(sorted(exclude_chromosomes)))
	#	sys.stdout.flush()
	print 'Regular expression used for filtering chromosome names:',chrom_pat.pattern

	# for statistics
	types = Counter()
	sources = Counter()

	# primary information
	genes = Counter()
	gene_chroms = dict()
	gene_ids = dict()

	# secondary information
	genes2 = Counter()
	polymorphic = set()

	# list of chromosomes
	chromosomes = set()

	i = 0
	missing = 0
	excluded_chromosomes = set()
	with open_plain_or_gzip(input_file) if input_file != '-' else sys.stdin as fh:
		#if i >= 500000: break
		reader = csv.reader(fh,dialect='excel-tab')
		for l in reader:
			i += 1
			if i % int(1e5) == 0:
				print '\r%d...' %(i), ; sys.stdout.flush() # report progress

			if len(l) > 1 and l[2] == field_name:

				# Note: Older Ensembl GTF files sometimes have a leading space in field #9
				attr = parse_attributes(l[8].lstrip(' '))

				type_ = attr['gene_biotype']
				if type_ not in ['protein_coding','polymorphic_pseudogene']:
					continue

				chrom = l[0]

				# test whether chromosome is valid
				m = chrom_pat.match(chrom)
				if m is None:
					excluded_chromosomes.add(chrom)
					continue

				chromosomes.add(m.group())

				source = l[1]
				id_ = attr['gene_id']
				try:
					name = attr['gene_name']
				except KeyError as e:
					missing += 1
					continue

				# store gene name
				genes[name] += 1

				# store Ensemble ID
				try:
					gene_ids[name].add(id_)
				except KeyError:
					gene_ids[name] = set([id_])

				# store chromsome
				try:
					gene_chroms[name].add(chrom)
				except KeyError:
					gene_chroms[name] = set([chrom])

				# record some statistics
				sources[source] += 1
				types[type_] += 1
				if type_ == 'polymorphic_pseudogene':
					polymorphic.add(name)
					genes2[name] += 1

	print "done (parsed %d lines)." %(i)

	print 
	print "Gene chromosomes (%d):" %(len(chromosomes))
	print "\t" + ', '.join(sorted(chromosomes))
	print
	print "Excluded chromosomes (%d):" %(len(excluded_chromosomes))
	print "\t" + ', '.join(sorted(excluded_chromosomes))

	if missing > 0:
		print 
		print "Genes without names:",missing

	def print_counter(C,sep=os.linesep):
		print sep.join(['\t%s: %d' %(k,C[k]) for k in sorted(C.keys(),key=lambda x:-C[x])])

	print
	print "Gene sources:"
	print_counter(sources)

	print
	print "Gene types:"
	print_counter(types)

	redundant_genes = sorted(g for g in genes if genes[g]>1)
	print
	print 'Genes with redundant annotations:', len(redundant_genes)

	print
	print 'Polymorphic pseudogenes (%d):' %(len(polymorphic)), ', '.join('%s (%d)' %(g,genes2[g]) for g in sorted(polymorphic))

	print
	print "Total protein-coding genes:", len(genes)

	with open(args.output_file,'w') as ofh:
		writer = csv.writer(ofh,dialect='excel-tab',lineterminator='\n',quoting=csv.QUOTE_NONE)
		for name in sorted(genes):
			chroms = ','.join(sorted(gene_chroms[name]))
			ids = ','.join(sorted(gene_ids[name]))
			writer.writerow([name,chroms,ids])

	return 0

if __name__ == '__main__':
	return_code = main()
	sys.exit(return_code)
