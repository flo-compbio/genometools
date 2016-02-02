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

"""Module containing the `ExpGenome` class."""

import os
import io
import logging
from collections import OrderedDict, Iterable

import unicodecsv as csv
import numpy as np

from genometools import misc
from .gene import ExpGene

logger = logging.getLogger(__name__)

class ExpGenome(object):
    """A complete set of genes in a gene expression analysis.

    The class represents a "genome" in the form of an ordered set of genes.
    This means that each gene has an index value, i.e. an integer indicating its
    0-based position in the genome.

    Parameters
    ----------
    genes: list or tuple of ExpGene objects
        The genes in the analysis.

    Notes
    -----
    The implementation is very similar to the `genometools.basic.GeneSetDB`
    class. It uses ordered dictionaries to support efficient access by gene
    name or index, as well as looking up the index of specific gene.
    """
    def __init__(self, genes):

        assert isinstance(genes, (list, tuple))
        for g in genes:
            assert isinstance(g, ExpGene)

        self._genes = OrderedDict([g.name, g] for g in genes)
        self._gene_names = tuple(self._genes.keys())
        self._gene_indices = OrderedDict([g.name, i] for i, g in enumerate(genes))
        logger.debug('Initialized ExpGenome with %d genes.', self.p)

    def __repr__(self):
        return '<%s object (%d genes; hash = %d)>' \
                %(self.__class__.__name__, self.p, hash(self))

    def __str__(self):
        return '<%s object (%d genes)>' %(self.__class__.__name__, self.p)

    def __getitem__(self, key):
        """Simple interface for accessing genes.

        Depending on whether key is an integer or not, look up a gene
        either by index, or by name.
        """
        if isinstance(key, int):
            return self.get_by_index(key)
        else:
            return self.get_by_name(key)

    def __hash__(self):
        return hash(self.genes)

    def __eq__(self, other):
        if self is other:
            return True
        elif type(self) != type(other):
            return False
        else:
            return repr(self) == repr(other)

    def __ne__(self, other):
        return not (self == other)

    def __contains__(self, gene):
        return gene in self._genes

    @property
    def p(self):
        """Returns the number of genes."""
        return len(self._genes)

    @property
    def genes(self):
        """Returns a tuple with all gene names."""
        return tuple(self._genes.values())

    @classmethod
    def from_gene_names(cls, genes):
        """Generate a genome from a list of gene names.

        Parameters
        ----------
        genes: list or tuple of (str, unicode)
            The list of gene names.
        
        Returns
        -------
        ExpGenome
            The genome.
        """
        assert isinstance(genes, (list, tuple))
        for g in genes:
            assert isinstance(g, (str, unicode))

        genome = cls([ExpGene(g) for g in genes])
        return genome

    def get_by_name(self, name):
        """Look up a gene by its name.

        Parameters
        ----------
        name: str or unicode
            The gene name.
        
        Returns
        -------
        `genometools.expression.ExpGene`
            The gene.
        """
        try:
            return self._genes[name]
        except KeyError:
            raise ValueError('No gene with name "%s"!' %(id_))

    def get_by_index(self, i):
        """Look up a gene by its index.

        Parameters
        ----------
        i: int
            The index.
        
        Returns
        -------
        `genometools.expression.ExpGene`
            The gene.
        """
        if i >= self.p:
            raise ValueError('Index %d out of bounds '
                    'for genome with %d genes.' %(i, self.p))

        return self._genes[self._gene_names[i]]

    def index(self, gene_name):
        """Returns the index of the gene with the given gene.

        The index is 0-based, so the first gene in the genome has the index 0,
        and the last one has index p-1.

        Parameters
        ----------
        gene_name: str or unicode
            The gene name.

        Returns
        -------
        int
            The gene index.
        """
        assert isinstance(gene_name, (str, unicode))

        try:
            return self._gene_indices[gene_name]
        except KeyError:
            raise ValueError('No gene with name "%s"!' %(gene_name))

    @classmethod
    def read_tsv(cls, path, enc = 'UTF-8'):
        """Read genes from tab-delimited text file.

        Parameters
        ----------
        path: str
            The path of the text file.
        enc: str, optional
            The file encoding. ("UTF-8")

        Returns
        -------
        None
        """
        genes = []
        with io.open(path, 'rb') as fh:
            reader = csv.reader(fh, dialect = 'excel-tab', encoding = enc)
            for l in reader:
                genes.append(ExpGene.from_list(l))
        return cls(genes)

    def write_tsv(self, path, enc = 'UTF-8'):
        """Write genes to tab-delimited text file in alphabetical order.

        Parameters
        ----------
        path: str
            The path of the output file.
        enc: str, optional
            The file encoding. ("UTF-8")

        Returns
        -------
        None
        """
        with io.open(path, 'wb') as ofh:
            writer = csv.writer(ofh, dialect = 'excel-tab', encoding = enc,
                    lineterminator = os.linesep, quoting = csv.QUOTE_NONE)
            for g in genes:
                writer.writerow(g.to_list())
        logger.info('Wrote %d genes to file "%s".', self.p, path)
