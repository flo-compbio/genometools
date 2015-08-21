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

	#genes_with_tss = []

	#transcript_info = {}
	#transcript_ensembl = {}

	transcript_gene = {}
	transcript_chrom = {}
	transcript_tss = {}

	i = 0
	ignored = 0
	problem = False
	with open_plain_or_gzip(input_file) if input_file != '-' else sys.stdin as fh:
		reader = csv.reader(fh,dialect='excel-tab')
		for l in reader:

			if i % int(1e5) == 0:
				print '\r%d...' %(i), ; sys.stdout.flush() # report progress
			i += 1
			#if i >= 1000:break

			if len(l) > 1 and l[2] == field_name:
				attr = parse_attributes(l[8].lstrip(' '))
				type_ = attr['gene_biotype']

				if type_ in ['protein_coding','polymorphic_pseudogene']:

					chrom = l[0]

					# test whether chromosome is valid
					m = chrom_pat.match(chrom)
					if m is None:
						excluded_chromosomes.add(chrom)
						continue

					source = l[1]

					try:
						gene_id = attr['gene_id']
					except KeyError as e:
						print attr_sep.split(l[8])
						print i,l
						print attr
						raise e

					gene_name = attr['gene_name']
					transcript_name = attr['transcript_name']

					if gene_id in ignore_ensembl:
						ignored += 1
						continue

					pos = 0
					if l[6] == '+': pos = int(l[3])-1
					else: pos = -(int(l[4])-1)

					if transcript_name in transcript_chrom:

						# another exon of this transcript has been seen before
						# make sure transcript information is consistent
						tchrom = transcript_chrom[transcript_name]
						tg = transcript_gene[transcript_name]
						tpos = transcript_tss[transcript_name]
						assert tg == gene_name
						assert tchrom == chrom
						assert sign(tpos) == sign(pos)

						# if exon is located more upstream than previous exons of this transcript,
						# remember its position
						transcript_tss[transcript_name] = min(tpos,pos)
						
					else:
						transcript_chrom[transcript_name] = chrom
						transcript_gene[transcript_name] = gene_name
						transcript_tss[transcript_name] = pos

	
	print "done (parsed %d lines)." %(i)

	if problem:
		print 'There were transcripts that had multiple annotations with different orientations.'
		print 'Please resolve these conflicts by specifying Ensembl IDs to be excluded.'
		return 1

	#counts = np.int64([len(transcript_ensembl[t]) for t in all_transcripts])
	#print 'Ensembl Gene ID counts per transcript name:'
	#print np.bincount(counts)

	#counts = np.int64([len(tss[t]) for t in all_transcripts])
	#print 'Transcript annotation counts per transcript name:'
	#print np.bincount(counts)

	# output
	all_transcripts = sorted(transcript_chrom.keys())
	with open(output_file,'w') as ofh:
		writer = csv.writer(ofh,dialect='excel-tab',lineterminator=os.linesep,quoting=csv.QUOTE_NONE)
		for t in all_transcripts:
			tss = [t, transcript_gene[t], transcript_chrom[t], str(transcript_tss[t])]
			writer.writerow(tss)

	return 0

if __name__ == '__main__':
	return_code = main()
	sys.exit(return_code)
