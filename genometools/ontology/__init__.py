"""Sub-package for working with Gene Ontology data."""

#from __future__ import (absolute_import, division,
#                        print_function)
#from builtins import *

from .term import GOTerm
from .ontology import GeneOntology
from .annotation import GOAnnotation
from .gaf import parse_gaf, get_goa_gene_sets
from .util import get_current_ontology_date, download_release

#__all__ = ['GOTerm', 'GeneOntology',
#           'GOAnnotation', 'parse_gaf', 'get_goa_gene_sets']