#!/usr/bin/env python

# Note: Some transcripts occur in the annotation of multiple (Ensembl) genes, sometimes with opposing strand orientations

import sys
import os
import gzip
import csv
import re
import gzip
import argparse
from math import copysign

import numpy as np

from genometools import misc
from genometools import genomic_features

def read_args_from_cmdline():
	parser = argparse.ArgumentParser(description='')

	parser.add_argument('-a','--annotation-file',default='-')

	parser.add_argument('-o','--output-file',required=True)

	parser.add_argument('-c','--chromosome-pattern',default=r'(?:\d\d?|MT|X|Y)$')
	parser.add_argument('-f','--field-name',default='exon')
	parser.add_argument('-i','--ignore-ensembl',default=[],nargs='+')

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

sign = lambda x:int(copysign(1.0,x))

def main(args=None):

	if args is None:
		args = read_args_from_cmdline()

	input_file = args.annotation_file
	output_file = args.output_file
	chrom_pat = re.compile(args.chromosome_pattern)
	field_name = args.field_name
	ignore_ensembl = set(args.ignore_ensembl)

	# list of chromosomes
	chromosomes = set()
	excluded_chromosomes = set()

	tss = {}
	genes_with_tss = []

	transcript_info = {}
	transcript_ensembl = {}

	transcript_gene = {}
	transcript_chrom = {}
	transcript_tss = {}

	i = 0
	conflicts = 0
	problem = False
	with open_plain_or_gzip(input_file) if input_file != '-' else sys.stdin as fh:
		#if i >= 500000: break
		reader = csv.reader(fh,dialect='excel-tab')
		for l in reader:

			if i % int(1e5) == 0:
				print '\r%d...' %(i), ; sys.stdout.flush() # report progress
			i += 1

			if len(l) > 1 and [2] == field_name:
				attr = parse_attributes(l[8])
				type_ = attr['gene_biotype']

				if type_ in ['protein_coding','polymorphic_pseudogene']:

					chrom = l[0]

					# test whether chromosome is valid
					m = chrom_pat.match(chrom)
					if m is None:
						excluded_chromosomes.add(chrom)
						continue

					source = l[1]

					gene_id = attr['gene_id']
					gene_name = attr['gene_name']
					transcript_name = attr['transcript_name']


					if gene_id in ignore_ensembl:
						continue

					pos = 0
					if l[6] == '+': pos = int(l[3])-1
					else: pos = -(int(l[4])-1)

					try:
						tg = transcript_gene[transcript_name]
						assert tg == gene_name
						tpos = transcript_tss[transcript_name]
						tchrom = transcript_chrom[transcript_name]
						assert tchrom == chrom
						assert sign(tpos) == sign(pos)
						transcript_tss[transcript_name] = min(tpos,pos)
					except KeyError:
						transcript_gene = gene_name
						transcript_tss = pos

					try:
						transcript_ensembl[transcript_name].add(gene_id)
					except KeyError:
						transcript_ensembl[transcript_name] = set([gene_id])

					continue
					# 
					info = (gene_name,chrom,int(copysign(1,pos)))
					if transcript_name in transcript_info:
						try:
							assert transcript_info[transcript_name] == info
						except AssertionError as e:
							problem = True
							print
							print 'Problem:',transcript_name
							print transcript_info[transcript_name]
							print info
					else:
						transcript_info[transcript_name] = info

					start_site = (gene_name,chrom,str(pos))
					try:
						tss[transcript_name].append(start_site)
					except KeyError:
						tss[transcript_name] = [start_site]

	print "done (parsed %d lines)." %(i)

	if problem:
		print 'There were transcripts that had multiple annotations with different orientations.'
		print 'Please resolve these conflicts by specifying Ensembl IDs to be excluded.'
		return 1

	#all_transcripts = sorted(transcript_ensembl.keys())

	#counts = np.int64([len(transcript_ensembl[t]) for t in all_transcripts])
	#print 'Ensembl Gene ID counts per transcript name:'
	#print np.bincount(counts)

	return 0

	#counts = np.int64([len(tss[t]) for t in all_transcripts])
	#print 'Transcript annotation counts per transcript name:'
	#print np.bincount(counts)

	c="""
	# for transcripts that have multiple annotations, select TSS to be the most upstream TSS
	for t in all_transcripts:
		start_sites = tss[t]
		if int(start_sites[0][2]) >= 0:
			start_sites = sorted(start_sites,key=lambda s:int(s[2]))
		else:
			start_sites = sorted(start_sites,key=lambda s:-int(s[2]))
		tss[t] = start_sites[0]
	"""

	# output
	with open(output_file,'w') as ofh:
		writer = csv.writer(ofh,dialect='excel-tab',lineterminator=os.linesep,quoting=csv.QUOTE_NONE)
		for t in all_transcripts:
			writer.writerow([t] + list(tss[t]))

	return 0

if __name__ == '__main__':
	return_code = main()
	sys.exit(return_code)
