from genometools import misc

from HTSeq import GenomicPosition as GP, GenomicInterval as GI

class GenomicFeature(object):
	def __init__(self):
		pass
	
	def __repr__(self):
		return "<GenomicFeature>"
	def __str__(self):
		return "<GenomicFeature>"

	def __hash__(self):
		return hash(repr(self))

	def __eq__(self,other):
		if self is other:
			return True
		elif type(self) != type(other):
			return False
		elif repr(self) == repr(other):
			return True
		else:
			return False

class LocationFeature(GenomicFeature):
	def __init__(self,chrom,pos,strand='.'):
		GenomicFeature.__init__(self)
		assert strand in ['.','-','+']
		self.loc = GP(chrom,pos,strand)

	@property
	def chrom(self):
		return self.loc.chrom

	@property
	def pos(self):
		return self.loc.pos

	@property
	def strand(self):
		return self.loc.strand

	def __repr__(self):
		return '<LocationFeature (loc:"%s")>' %(repr(self.loc))

	def __str__(self):
		return '<LocationFeature at "%s">' %(str(self.loc))

	def set_location(self,chrom,pos,strand='.'):
		self.loc = GP(chrom,pos,strand)


class IntervalFeature(GenomicFeature):
	def __init__(self,chrom,start,end,strand='.'):
		assert start < end
		assert strand in ['.','-','+']
		self.iv = GI(chrom,start,end,strand)

	@property
	def chrom(self):
		return self.iv.chrom

	@property
	def start(self):
		return self.iv.start

	@property
	def end(self):
		return self.iv.end

	@property
	def strand(self):
		return self.iv.strand

	@property
	def length(self):
		return self.iv.length

	def __repr__(self):
		return '<IntervalFeature (iv:"%s")>' %(repr(self.iv))

	def __str__(self):
		return '<IntervalFeature with interval "%s">' %(str(self.iv))

	def set_interval(self,chrom,start,end,strand='.'):
		self.iv = GI(chrom,start,end,strand)


class GeneFeature(GenomicFeature):
	def __init__(self,gene):
		self.gene = gene

	def __repr__(self):
		return '<GeneFeature (gene:"%s")>' %(self.gene)

	def __str__(self):
		return '<GeneFeature of gene "%s">' %(self.gene)


class TranscriptFeature(GeneFeature):
	def __init__(self,gene,transcript):
		GeneFeature.__init__(self,gene)
		self.transcript = transcript

	def __repr__(self):
		return '<TranscriptFeature (gene:"%s", transcript:"%s")>' %(self.gene,self.transcript)

	def __str__(self):
		return '<TranscriptFeature of gene "%s", transcript "%s">' %(self.gene,self.transcript)


class ChipPeak(IntervalFeature):
	def __init__(self,chrom,start,end,summit):
		IntervalFeature.__init__(self,chrom,start,end)
		assert summit < (end-start)
		self.summit = summit

	@classmethod
	def read_ucsc_narrowpeaks(cls,narrowpeak_file):
		data = misc.read_all(narrowpeak_file)
		peaks = []
		for d in data:

			# convert chromosome name to Ensembl
			chrom = d[0]
			if chrom[:3] == 'chr':
				chrom = chrom[3:]
			if chrom == 'M':
				chrom = 'MT'

			P = cls(d[0],int(d[1]),int(d[2]),int(d[9]))
			peaks.append(P)

		return peaks

	def __repr__(self):
		return "<ChipPeak (chrom:%s,start:%d,end:%d,summit:%d)>" %(self.chrom,self.start,self.end,self.summit)

	def __str__(self):
		return "<ChipPeak on chromosome '%s' (%d - %d), length = %d bp, summit @ %d>" \
				%(self.chrom,self.start,self.end,self.length,self.summit+1)


class TSS(TranscriptFeature,LocationFeature):
	def __init__(self,gene,transcript,chrom,pos,strand):
		TranscriptFeature.__init__(self,gene,transcript)
		LocationFeature.__init__(self,chrom,pos,strand)

	def __repr__(self):
		return '<TSS (gene:"%s", transcript:"%s", strand:%d)>' %(self.gene,self.transcript,self.strand)
	def __str__(self):
		return '<TSS of gene "%s"/transcript "%s", on "%s" strand>' %(self.gene,self.transcript,self.strand)


def convert_chrom_from_ensembl(c):
	chrom = c
	try:
		chrom = int(chrom)
		chrom = 'chr' + str(chrom)
	except ValueError:
		if chrom == 'X':
			chrom = 'chrX'
		elif chrom == 'Y':
			chrom = 'chrY'
		elif chrom == 'MT':
			chrom = 'chrM'
	return chrom

def convert_chrom_to_ensembl(chrom):
	if chrom == 'chrM':
		chrom = 'MT'
	elif chrom.startswith('chr'):
		chrom = chrom[3:]
	return chrom
