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

import csv
import logging
from collections import OrderedDict

import numpy as np

from genometools import misc
from .gene import ExpGene

logger = logging.getLogger(__name__)

class ExpGenome(object):
    """A complete set of genes in a gene expression analysis.

    Parameters
    ----------
    genes: list or tuple of ExpGene objects
        The genes in the analysis.
    """
    def __init__(self, genes):
        genes = OrderedDict([g.name,g] for g in genes)

        # make sure genes are in alphabetical order
        names = genes.keys()
        a = np.lexsort([names])
        if not np.all(a == np.arange(len(genes))):
            sorted_names = [names[i] for i in a]
            genes = OrderedDict([n,genes[n]] for n in sorted_names)

        self._genes = genes
        logger.debug('Initialized ExpGenome with %d genes.', self.p)

    def __repr__(self):
        return '<ExpGenome (%d genes; hash = %d)>' \
                %(self.p, hash(self.get_genes()))

    def __str__(self):
        return '<ExpGenome with %d genes>' %(self.p)

    def __hash__(self):
        return hash(repr(self))

    def __eq__(self,other):
        if type(self) != type(other):
            return False
        elif repr(self) == repr(other):
            return True
        else:
            return False

    @property
    def p(self):
        return len(self._genes)

    def get_genes(self):
        return tuple(self._genes.values())

    def get_gene_names(self):
        return self._genes.keys()

    def index(self,gene):
        return misc.bisect_index(self.get_gene_names(),gene)

    def has_gene(self,name):
        return name in self._genes

    def get_gene(self,name):
        return self._genes[name]

    @classmethod
    def read_tsv(cls,path):
        """Read genes from tab-delimited text file.

        Parameters
        ----------
        path: str
            The path of the text file.
        """
        genes = []
        with open(path) as fh:
            reader = csv.reader(fh, dialect = 'excel-tab')
            for l in reader:
                genes.append(ExpGene.from_list(l))
        return cls(genes)

    def write_tsv(self,path):
        """Write genes to tab-delimited text file in alphabetical order.

        Parameters
        ----------
        path: str
            The path of the output file.
        """
        with open(path,'w') as ofh:
            writer = csv.writer(ofh, dialect = 'excel-tab',
                    lineterminator = '\n', quoting = csv.QUOTE_NONE)
            for g in self._genes.itervalues():
                writer.writerow(g.to_list())
        logger.info('Wrote %d genes to file "%s".', self.p, path)
