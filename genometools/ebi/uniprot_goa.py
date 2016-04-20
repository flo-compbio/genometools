# Copyright (c) 2015, 2016 Florian Wagner
#
# This file is part of GenomeTools.
#
# GenomeTools is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License, Version 3,
# as published by the Free Software Foundation.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

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