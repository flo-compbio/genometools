import pkg_resources

__version__ = pkg_resources.require('genometools')[0].version



import genometools.misc

from genometools import gtf, ensembl, rnaseq, sra

__all__ = []
