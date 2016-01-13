# Copyright (c) 2015 Florian Wagner
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

"""Module containing the `ExpGenome` class."""

import os
import logging
from collections import OrderedDict, Iterable

import unicodecsv as csv
import numpy as np

from genometools import misc
from .gene import ExpGene

logger = logging.getLogger(__name__)

class ExpGenome(object):
    """A complete set of genes in a gene expression analysis.

    Parameters
    ----------
    genes: list (tuple, set) of ExpGene objects
        The genes in the analysis.
    """
    def __init__(self, genes):

        assert isinstance(genes, Iterable)

        genes = sorted(genes, key = lambda g: g.name)
        self._genes = OrderedDict([g.name, g] for g in genes)
        self._gene_indices = OrderedDict([g.name, i] for i, g in enumerate(genes))
        logger.debug('Initialized ExpGenome with %d genes.', self.p)

    def __repr__(self):
        return '<ExpGenome (%d genes; hash = %d)>' \
                %(self.p, hash(self))

    def __str__(self):
        return '<ExpGenome with %d genes>' %(self.p)

    def __hash__(self):
        return hash(self.genes)

    def __eq__(self,other):
        if self is other:
            return True
        elif type(self) != type(other):
            return False
        else:
            return  repr(self) == repr(other):

    def __contains__(self, gene):
        return gene in self._genes

    @property
    def p(self):
        """Returns the number of genes."""
        return len(self._genes)

    def genes(self):
        """Returns the list of genes as a tuple."""
        return tuple(self._genes.values())

    def index(self, gene):
        """Returns the index of the given gene."""
        return self.gene_indices[gene]

    @classmethod
    def read_tsv(cls, path, enc = 'utf-8'):
        """Read genes from tab-delimited text file.

        Parameters
        ----------
        path: str
            The path of the text file.
        enc: str
            The file encoding.
        """
        genes = []
        with open(path, 'rb') as fh:
            reader = csv.reader(fh, dialect = 'excel-tab', encoding = enc)
            for l in reader:
                genes.append(ExpGene.from_list(l))
        return cls(genes)

    def write_tsv(self, path, enc = 'utf-8'):
        """Write genes to tab-delimited text file in alphabetical order.

        Parameters
        ----------
        path: str
            The path of the output file.
        enc: str
            The file encoding.
        """
        with open(path, 'wb') as ofh:
            writer = csv.writer(ofh, dialect = 'excel-tab', encoding = enc,
                    lineterminator = os.linesep, quoting = csv.QUOTE_NONE)
            for g in self._genes.itervalues():
                writer.writerow(g.to_list())
        logger.info('Wrote %d genes to file "%s".', self.p, path)
