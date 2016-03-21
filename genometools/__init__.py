import os
import pkg_resources

__version__ = pkg_resources.require('genometools')[0].version

_root = os.path.abspath(os.path.dirname(__file__))

from genometools import misc, gtf, ensembl, rnaseq, sra

__all__ = [misc, gtf, ensembl, rnaseq, rnaseq]
