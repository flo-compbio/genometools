import pkg_resources

__version__ = pkg_resources.require('genometools')[0].version

from genometools import misc, gtf, ensembl, rnaseq, sra

__all__ = []
