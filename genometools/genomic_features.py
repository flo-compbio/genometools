import sys
import os
import csv
import gzip
from math import copysign

from genometools import misc

from HTSeq import GenomicPosition as GP, GenomicInterval as GI

sign = lambda x:int(copysign(1.0,x))

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

	def get_interval(self):
		return GI(self.chrom,self.pos,self.pos+1,self.strand)


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
	def __init__(self,transcript,gene):
		GeneFeature.__init__(self,gene)
		self.transcript = transcript

	def __repr__(self):
		return '<TranscriptFeature (transcript:"%s", gene:"%s")>' %(self.transcript,self.gene)

	def __str__(self):
		return '<TranscriptFeature of transcript "%s", gene "%s">' %(self.transcript,self.gene)


class ChipPeak(IntervalFeature):
	def __init__(self,exp,chrom,start,end,summit=None):
		IntervalFeature.__init__(self,chrom,start,end)
		if summit is not None:
			assert summit >= 0
			assert summit < (end-start)
		self.summit = summit
		self.exp = exp

	def __repr__(self):
		return '<ChipPeak (exp:%s, chrom:%s, start:%d, end:%d, summit:%s)>' \
				%(self.exp,self.chrom,self.start,self.end,str(self.summit))

	def __str__(self):
		summit_str = 'summit N/A'
		if self.summit is not None:
			summit_str = 'summit @ %d' %(self.summit)
		return '<ChipPeak from experiment "%s" on chromosome "%s" (%d - %d), length = %d bp, %s>' \
				%(self.exp,self.chrom,self.start,self.end,self.length,summit_str)

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

	def read_ucsc_bed_regions(cls,file_name,experiment=None):
		# reads peaks from bed file (no summit)
		data = misc.read_all(file_name)
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

	def get_interval(self):
		return GI(self.chrom,self.start,self.end)

	def get_summit_interval(self):
		if self.summit is None:
			raise ValueError('This peak does not have a summit position specified.')

		return GI(self.chrom,self.start+self.summit,self.start+self.summit+1)


class GeneIntervalFeature(GeneFeature,IntervalFeature):
	def __init__(self,gene,chrom,start,end,strand='.'):
		GeneFeature.__init__(self,gene)
		IntervalFeature.__init__(self,chrom,start,end,strand)

	def __repr__(self):
		return '<GeneIntervalFeature (gene:"%s", chrom:"%s", strand:"%s", loc:[%d;%d))>' %(self.gene,self.chrom,self.strand,self.start,self.end)

	def __str__(self):
		return '<GeneIntervalFeature of gene "%s" on chromosome "%s", strand "%s", location [%d;%d)>' %(self.gene,self.chrom,self.strand,self.start,self.end)

class TranscriptIntervalFeature(TranscriptFeature,IntervalFeature):
	def __init__(self,transcript,gene,chrom,start,end,strand='.'):
		TranscriptFeature.__init__(self,transcript,gene)
		IntervalFeature.__init__(self,chrom,start,end,strand)

	def __repr__(self):
		return '<TranscriptIntervalFeature (transcript: "%s", gene:"%s", chrom:"%s", strand:"%s", loc:[%d;%d))>' \
				%(self.transcript,self.gene,self.chrom,self.strand,self.start,self.end)

	def __str__(self):
		return '<TranscriptIntervalFeature of transcript "%s"/gene "%s" on chromosome "%s", strand "%s", location [%d;%d)>' \
				%(self.transcript,self.gene,self.chrom,self.strand,self.start,self.end)


class Promoter(GenomicFeature):
	""" 
	Special Feature that contains a TSS object (with strand information),
	as well as the distance upstream and downstream of it (in bp) that
	should be considered the promoter.
	"""

	def __init__(self,tss,upstream,downstream):
		GenomicFeature.__init__(self)
		assert tss.strand in ['+','-']
		self.tss = tss
		self.upstream = upstream
		self.downstream = downstream

	@property
	def transcript(self):
		return self.tss.transcript

	@property
	def gene(self):
		return self.tss.gene

	@property
	def chrom(self):
		return self.tss.chrom

	@property
	def pos(self):
		return self.tss.pos

	@property
	def strand(self):
		return self.tss.strand

	def __repr__(self):
		return '<Promoter (tss:"%s", upstream:%d, downstream:%d>' \
				%(repr(self.tss),self.upstream,self.downstream)

	def __str__(self):
		return '<Promoter of %s, [-%d;%d)>' %(str(self.tss),self.upstream,self.downstream)

	def get_interval(self,chromlen=None,stranded=False):
		""" Returns the full interval, containing upstream and downstream portions of the promoter. """
		pos = self.pos

		left = None
		right = None
		if self.strand == '+':
			left = pos - self.upstream
			right = pos + self.downstream
		else:
			right = pos + self.upstream + 1
			left = pos - self.downstream + 1
		if chromlen is not None:
			assert self.chrom in chromlen
			left = max(0,left)
			right = min(right,chromlen[self.chrom])
		iv = None
		strand = self.strand
		if not stranded:
			strand = '.'
		iv = GI(self.chrom,left,right,strand)
		return iv


class TSS(TranscriptFeature,LocationFeature):

	def __init__(self,transcript,gene,chrom,pos,strand):
		TranscriptFeature.__init__(self,transcript,gene)
		LocationFeature.__init__(self,chrom,pos,strand)

	orient = {1: '+', -1: '-'}

	def __repr__(self):
		return '<TSS (transcript:"%s", gene:"%s", loc:%s)>' %(self.transcript,self.gene,repr(self.loc))

	def __str__(self):
		return '<TSS of transcript "%s"/gene "%s", at "%s">' %(self.transcript,self.gene,str(self.loc))

	@classmethod
	def read_tss_file(cls,fn):
		orient = cls.orient
		start_sites = []
		with open(fn) as fh:
			reader = csv.reader(fh,dialect='excel-tab')
			for l in reader:
				o = orient[sign(int(l[3]))]
				pos = abs(int(l[3]))
				start_sites.append(cls(l[0],l[1],l[2],pos,o))
		return start_sites
		

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
