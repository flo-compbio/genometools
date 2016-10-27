"""Sub-package for working with Gene Ontology data."""

#from __future__ import (absolute_import, division,
#                        print_function)
#from builtins import *

from .term import GOTerm
from .ontology import GeneOntology
from .annotation import GOAnnotation
from .gaf import parse_gaf 

__all__ = ['GOTerm', 'GeneOntology',
           'GOAnnotation', 'parse_gaf']