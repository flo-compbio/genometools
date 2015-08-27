# generate tss data structures
import sys
import re

import numpy as np

from HTSeq import GenomicArrayOfSets, GenomicInterval as GI

from genometools.features import IntervalFeature,TSS
from genometools import misc

#chromlen_file = output_dir + 'chromosome_lengths_mouse_legacy.tsv'
#tss_file = output_dir + 'protein_coding_tss_mouse_legacy.tsv'

#chrom_pat = re.compile(r'(?:\d\d?|X)$')

class ChipPeak(IntervalFeature):
	def __init__(self,exp,chrom,start,end,summit=None,pval=None):
		IntervalFeature.__init__(self,chrom,start,end)
		if summit is not None:
			assert summit >= 0
			assert summit < (end-start)
		if pval is not None:
			assert pval >= 0
		self.summit = summit
		self.exp = exp
		self.pval = pval

	def __repr__(self):
		return '<ChipPeak (exp:%s, chrom:%s, start:%d, end:%d, summit:%s, pval:%s)>' \
				%(self.exp,self.chrom,self.start,self.end,str(self.summit),str(self.pval))

	def __str__(self):
		summit_str = 'summit N/A'
		if self.summit is not None:
			summit_str = 'summit @ %d' %(self.summit)
		pval_str = 'pval N/A'
		if self.pval is not None:
			pval_str = 'pval = %.1e' %(self.pval)
		return '<ChipPeak from experiment "%s" on chromosome "%s" (%d - %d), length = %d bp, %s, %s>' %(summit_str,pval_str) \
				%(self.exp,self.chrom,self.start,self.end,self.length,summit_str,self.pval)

	@staticmethod
	def get_chrom_ucsc2ensembl(chrom):
		if chrom[:3] == 'chr':
			chrom = chrom[3:]
		if chrom == 'M':
			chrom = 'MT'
		return chrom

	@staticmethod
	def get_chrom_ensembl2ucsc(chrom):
		if chrom == 'MT':
			chrom = 'M'
		chrom = 'chr' + chrom
		return chrom

	@classmethod
	def from_text(cls,text,sep='\t'):
		data = text.split(sep)
		summit = None
		if data[3] != 'None':
			summit = int(data[3])
		return cls(data[0],int(data[1]),int(data[2]),summit)

	@classmethod
	def read_encode_narrowpeaks(cls,narrowpeak_file,experiment=None):
		# assumes that summit positions are specified (last column)
		data = misc.read_all(narrowpeak_file)
		peaks = []
		for d in data:
			# convert chromosome name to Ensembl
			chrom = ChipPeak.get_chrom_ucsc2ensembl(d[0])
			try:
				summit = int(d[9])-1 # -1 since Encode apparently isn't using format correctly
				P = cls(experiment,chrom,int(d[1]),int(d[2]),summit)
			except AssertionError as e:
				print chrom,d[1],d[2],d[9]
				raise e
			peaks.append(P)

		return peaks

	@classmethod
	def read_ucsc_narrowpeaks(cls,narrowpeak_file,experiment=None):
		# assumes that summit positions are specified (last column)
		data = misc.read_all(narrowpeak_file)
		peaks = []
		for d in data:
			# convert chromosome name to Ensembl
			chrom = ChipPeak.get_chrom_ucsc2ensembl(d[0])
			pval = pow(10,-float(d[-3])/10.0)
			try:
				summit = int(d[9])-1 # -1 since Encode apparently isn't using format correctly
				P = cls(experiment,chrom,int(d[1]),int(d[2]),summit,pval)
			except AssertionError as e:
				print chrom,d[1],d[2],d[9],d[-3]
				raise e
			peaks.append(P)

		return peaks

	@classmethod
	def read_encode_narrowpeaks_fixed_width(cls,narrowpeak_file,chromlen,width,experiment=None):
		# assumes that summit positions are specified (last column)
		data = misc.read_all(narrowpeak_file)
		goleft = int(width/2.0)
		goright = width-goleft
		peaks = []
		for d in data:
			chrom = ChipPeak.get_chrom_ucsc2ensembl(d[0])
			assert chrom in chromlen
			start = int(d[1])

			summit = int(d[9])-1
			pos = start + summit # genomic position of the peak summit

			left = max(pos-goleft,0)
			right = min(pos+goright,chromlen[chrom])
			P = cls(experiment,chrom,left,right,pos-left)
			peaks.append(P)

		return peaks

	@classmethod
	def read_ucsc_bed_regions(cls,file_name,experiment=None):
		# reads peaks from bed file (no summit)
		data = misc.read_all(file_name)
		data = data[1:] # remove header
		peaks = []
		for d in data:
			# convert chromosome name to Ensembl
			chrom = ChipPeak.get_chrom_ucsc2ensembl(d[0])
			P = cls(experiment,chrom,int(d[1]),int(d[2]))
			peaks.append(P)
		return peaks

	@staticmethod
	def write_encode_narrowpeaks(peaks,output_file,compressed=False):
		# writes to encode narrowpeaks format, makes summit indexing 1-based (violating the format)
		with open(output_file,'w') if not compressed else gzip.open(output_file,'w') as ofh:
			writer = csv.writer(ofh,dialect='excel-tab',lineterminator=os.linesep,quoting=csv.QUOTE_NONE)
			for p in peaks:
				chrom = ChipPeak.get_chrom_ensembl2ucsc(p.chrom)
				summit_str = str(-1)
				if p.summit is not None:
					summit_str = str(p.summit+1)
				writer.writerow([chrom,str(p.start),str(p.end),'.',str(0),'.',str(0),str(-1),str(-1),summit_str])

	def get_text(self,sep='\t'):
		return sep.join([self.chrom,str(self.start),str(self.end),str(self.summit)])

	def get_summit_position(self):
		if self.summit is None:
			raise ValueError('This peak does not have a summit position specified.')

		pos = self.iv.start + self.summit
		return GI(self.chrom,pos,pos+1)

	def get_summit_location(self):
		if self.summit is None:
			raise ValueError('This peak does not have a summit position specified.')

		return self.iv.start + self.summit

	def get_interval(self):
		return GI(self.chrom,self.start,self.end)

	def get_summit_interval(self):
		if self.summit is None:
			raise ValueError('This peak does not have a summit position specified.')

		return GI(self.chrom,self.start+self.summit,self.start+self.summit+1)


def get_genome_tss_cover(chromlen_file,tss_file,chrom_pat=None):

	# in order to efficiently calculate peak-gene distances, we need a "TSS cover"

	# read chromosome lengths
	chromlen = misc.read_chromlen(chromlen_file)

	# read transcription start sites
	start_sites = TSS.read_tss_file(tss_file)

	# filtering by chromosome
	if chrom_pat is not None:
		chromlen = dict([chrom,l] for chrom,l in chromlen.iteritems() if re.match(chrom_pat,chrom))

		# filter TSS
		start_sites = [s for s in start_sites if re.match(chrom_pat,s.chrom) is not None]

	print 'Annotating %d chromosomes with %d transcription start sites...' \
			%(len(chromlen),len(start_sites))
	sys.stdout.flush()

	genome_tss = GenomicArrayOfSets(chromlen,stranded=False)

	print 'Step 1/2: Annotating genome with TSS...',
	sys.stdout.flush()
	for s in start_sites:
		if s.chrom in chromlen:
			genome_tss[s.get_interval()] += s
	print 'done!'

	print 'Step 2/2: Generating genome cover...',
	sys.stdout.flush()
	genome_tss_cover = GenomicArrayOfSets(chromlen,stranded=False)
	for chrom in sorted(chromlen.keys()):
		sys.stdout.flush()
		prev_features = set()
		prev_iv = None
		for iv,features in genome_tss.chrom_vectors[chrom]['.'].steps():
			# go over all transcript start sites
			if features:
				# we hit one (or more) TSS
				if prev_iv is not None:
					# add this/these TSS to the previous inter-TSS interval
					for f in list(features):
						genome_tss_cover[prev_iv] += f
				assert iv.length == 1
				# at the current position, the current TSS is necessarily "the closest"
				for f in list(features):
					genome_tss_cover[iv] += f
				# store this TSS for future reference
				prev_features = features
			else:
				# we're in an inter-TSS interval
				if prev_features:
					# add the previous TSS
					for f in list(prev_features):
						genome_tss_cover[iv] += f
				# store this interval for future reference
				prev_iv = iv
		# the last interval was the one following the last TSS of the chromosome
		# when we encountered that interval, we added the previous TSS
		# therefore, we're done
	print 'done!'
	sys.stdout.flush()

	return genome_tss_cover


def get_genome_peak_cover(chromlen_file,peak_file):

	# in order to efficiently calculate gene-peak distances, we need a "peak cover"
	# we assume that peaks are non-overlapping

	from genometools.features import ChipPeak
	import common

	#chromlen_file = output_dir + 'chromosome_lengths_mouse.tsv'
	#demo_dir = data_dir + 'great_demo' + os.sep
	#peak_file = demo_dir + 'limb.mm10.bed'

	#chrom_pat = re.compile(r'(?:\d\d?|X)')

	# read chromosome lengths
	chromlen = common.read_chromlen(chromlen_file)

	# read peaks
	peaks = ChipPeak.read_ucsc_bed_regions(peak_file)

	# make a GenomicArray that holds all the peaks
	print 'Storing peaks in GenomicArray...',
	sys.stdout.flush()
	genome_peaks = GenomicArrayOfSets(chromlen,stranded=False)
	for p in peaks:
		if p.chrom in chromlen:
			genome_peaks[p.get_interval()] += p
	print 'done!'
	sys.stdout.flush()

	# use the peak GenomicArray to find all inter-peak intervals and annotate them appropriately
	print 'Generating peak cover...',
	sys.stdout.flush()
	genome_peak_cover = GenomicArrayOfSets(chromlen,stranded=False)
	for chrom in sorted(chromlen.keys()):
		prev_peak = None
		prev_iv = None
		for iv,features in genome_peaks.chrom_vectors[chrom]['.'].steps():
			# go over all peaks on this chromosome
			if features:
				# we hit a peak
				assert len(features) == 1 # we assume peaks are non-overlapping
				p = list(features)[0]
				if prev_iv is not None:
					# add this peak to the previous inter-peak interval
					# the previous inter-peak interval should now be annotated with its two adjacent peaks
					genome_peak_cover[prev_iv] += p
				# at the current position, the current peak is necessarily "the closest"
				genome_peak_cover[iv] += p
				# store this peak for future reference
				prev_peak = p
			else:
				# we're in an inter-peak interval
				if prev_peak is not None:
					# add the previous peak
					genome_peak_cover[iv] += prev_peak
				# store this interval for future reference
				prev_iv = iv
		# the last interval was the one following the last peak of the chromosome
		# when we encountered that interval, we added the previous peak
		# therefore, we're done
	print 'done!'
	sys.stdout.flush()

	return genome_peak_cover


def get_peak_gene_dists(genome_tss_cover,peaks,use_summit=False):
	all_chromosomes = set(genome_tss_cover.chrom_vectors.keys())
	n = len(peaks)
	peak_dists = np.zeros(n,dtype=np.int64)
	peak_dists += int(1e9)
	ignored = 0
	ignored_chroms = set()
	for i,p in enumerate(peaks):
		if p.chrom not in all_chromosomes:
			ignored += 1
			ignored_chroms.add(p.chrom)
			peak_dists[i] = float('nan')
			continue

		peak_iv = None
		if use_summit:
			peak_iv = p.get_summit_interval()
		else:
			peak_iv = p.get_interval()

		dist = int(1e9)
		for iv,features in genome_tss_cover[peak_iv].steps():
			# features are TSS objects
			for f in features:
				if f.pos >= peak_iv.start and f.pos < peak_iv.end:
					# TSS is contained within peak region
					# therefore the gene that this TSS belongs to has distance 0
					dist = 0
				else:
					# TSS is not contained within peak region
					# get closest distance from both peak boundaries
					dist = min(dist,min(abs(f.pos-peak_iv.start),abs(f.pos-(peak_iv.end-1))))
		peak_dists[i] = dist

	if ignored_chroms:
		print 'Warning: Ignored %d peaks on %d unknown chromosomes:' %(ignored,len(ignored_chroms))
		print ', '.join(sorted(ignored_chroms))

	return peak_dists

def get_gene_peak_dists(genome_peak_cover,genes,start_sites,use_summit=False):
	# we assume genes are in alphabetical order
	# go over all TSS and find minimum gene-peak distances
	print 'Determining minimum gene-peak distances...',
	sys.stdout.flush()
	n = len(genes)
	gene_dists = np.empty(n,dtype=np.int64)
	gene_dists[:] = int(1e9)
	m = len(start_sites)
	for i,s in enumerate(start_sites):
		try:
			idx = misc.bisect_index(genes,s.gene)
		except ValueError:
			continue
		iv = s.get_interval()
		dist = int(1e9)
		for iv,features in genome_peak_cover[iv].steps(): # features = peaks
			assert len(features) > 0 # make sure the whole genome is covered
			for p in features:
				if use_summit:
					# distance calculation is simple
					dist = min(dist,abs(s.pos-p.get_summit.location()))
				else:
					# two cases
					if s.pos >= p.start and s.pos < p.end:
						# TSS is contained within peak region
						# therefore the gene that this TSS belongs to has peak distance 0
						dist = 0
					else:
						# TSS is not contained within peak region
						# use minimum distance from both peak edges
						dist = min(dist,min(abs(s.pos-p.start),abs(s.pos-(p.end-1))))
		gene_dists[idx] = min(gene_dists[idx],dist)
	print 'done!'
	sys.stdout.flush()

	return gene_dists
