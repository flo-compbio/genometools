from genometools import misc

class GenomicFeature(object):
	def __init__(self):
		pass
	
	def __repr__(self):
		return "<GenomicFeature>" %(self.gene)
	def __str__(self):
		return "<GenomicFeature>" %(self.gene)

class ChipPeak(GenomicFeature):
	def __init__(self,chrom,start,end,summit):
		super(ChipPeak,self).__init__()
		assert start < end
		assert summit < (end-start)
		self.chrom = chrom
		self.start = start
		self.end = end
		self.summit = summit

	@property
	def length(self):
		return self.end - self.start

	@classmethod
	def read_ucsc_narrowpeak_file(cls,narrowpeak_file):
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
		return "<Peak %s: %d-%d (%d bp), summit=%d)>" %(self.chrom,self.start,self.end,self.length,self.summit)
	def __str__(self):
		return "<Peak on chromosome '%s' (%d - %d), length = %d bp, summit @ %d>" \
				%(self.chrom,self.start,self.end,self.length,self.summit+1)
		

class Peak(GenomicFeature):

	def __init__(self,chromosome,start,end,summit):
		super(Peak,self).__init__()
		assert start < end
		assert summit < (end-start)
		self.chromosome = chromosome
		self.start = start
		self.end = end
		self.summit = summit

	def get_length(self):
		return self.end - self.start

	def __repr__(self):
		return "<Peak %s: %d-%d (%d)>" %(self.chromosome,self.start,self.end,self.summit)
	def __str__(self):
		return "<Peak on chromosome '%d' (%d - %d), length = %d bp, summit @%d>" \
				%(self.chromosome,self.start,self.end,self.get_length(),self.summit+1)

	def __eq__(self,other):
		if type(self) != type(other):
			return False

		elif self is other:
			return True

		elif repr(self) == repr(other):
			return True

		else:
			return False
		"""			
		elif self.chromosome == other.chromosome \
			and self.start == other.start \
			and self.end == other.end \
			and self.summit == other.summit:
				return True

		else:
			return False
		"""

	def __hash__(self):
		return hash(repr(self))

"""
class MACSPeak(GenomicFeature):
	def __init__(self,data):
		super(MACSPeak,self).__init__()
		self.data = tuple(data) # tuple with coordinates etc., just like in MACS peak file

	def __repr__(self):
		return "<Feature gene:'%s'>" %(self.gene)
	def __str__(self):
		return "<Feature of gene '%s'>" %(self.gene)

	def __eq__(self,other):
		if isinstance(other,self.__class__):
			return (self.data == other.data)
		else:
			return False
	def __ne__(self,other):
		return not self.__eq__(other)

	def __hash__(self):
		return hash(self.data)
"""

class GeneFeature(GenomicFeature):
	def __init__(self,gene):
		super(GeneFeature,self).__init__()
		self.gene = gene

	def __repr__(self):
		return "<Feature gene:'%s'>" %(self.gene)
	def __str__(self):
		return "<Feature of gene '%s'>" %(self.gene)

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

class DirectedGeneFeature(GeneFeature):
	def __init__(self,gene,direction):
		super(DirectedGeneFeature,self).__init__(gene)
		self.dir = direction

	def __repr__(self):
		return "<DirectedGeneFeature (gene:'%s', direction %d)>" %(self.gene,self.dir)
	def __str__(self):
		return "<DirectedGeneFeature of gene '%s', direction %d>" %(self.gene,self.dir)

	def __eq__(self,other):
		if isinstance(other,self.__class__):
			return (self.gene == other.gene and self.dir == other.dir)
		else:
			return False
	def __ne__(self,other):
		return not self.__eq__(other)

	def __hash__(self):
		return hash((self.gene,self.dir))

class TSS(DirectedGeneFeature):
	def __init__(self,gene,transcript,direction):
		super(TSS,self).__init__(gene,direction)
		self.transcript = transcript

	def __repr__(self):
		return "<TSS (gene:'%s', transcript:'%s', direction %d)>" %(self.gene,self.transcript,self.dir)
	def __str__(self):
		return "<TSS of gene '%s'/transcript '%s', direction %d>" %(self.gene,self.transcript,self.dir)

	def __eq__(self,other):
		if not isinstance(other,self.__class__):
			return False
		elif other is self:
			return True
		elif repr(self) == repr(other):
			return True
			#return (self.gene == other.gene and self.transcript == other.transcript and self.dir == other.dir)
		else:
			return False
	def __ne__(self,other):
		return not self.__eq__(other)

	def __hash__(self):
		return hash((self.gene,self.transcript,self.dir))

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
