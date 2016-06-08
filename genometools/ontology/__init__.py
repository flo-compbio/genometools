"""Package for working with Gene Ontology data."""

from __future__ import (absolute_import, division,
                        print_function)
from builtins import *

from .term import GOTerm
from .ontology import GeneOntology

__all__ = ['GOTerm', 'GeneOntology']