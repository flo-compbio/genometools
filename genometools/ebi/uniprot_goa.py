"""Functions for working with UniProt-GOA data."""

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import *

# import io

from genometools import misc

def get_gaf_gene_ontology_file(path):
    """Extract the gene ontology file associated with a GO annotation file.

    Parameters
    ----------
    path: str
        The path name of the GO annotation file.

    Returns
    -------
    str
        The URL of the associated gene ontology file.
    """
    assert isinstance(path, str)

    version = None
    with misc.smart_open_read(path, encoding='UTF-8', try_gzip=True) as fh:
        for l in fh:
            if l[0] != '!':
                break
            if l.startswith('!GO-version:'):
                version = l.split(' ')[1]
                break
    return version