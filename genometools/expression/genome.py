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

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
_oldstr = str
from builtins import *

import os
import logging
import hashlib
from collections import OrderedDict, Iterable, Counter
import copy

import six
import pandas as pd
import unicodecsv as csv
import numpy as np

from .gene import ExpGene

logger = logging.getLogger(__name__)


class ExpGenome(object):
    """A complete set of genes in a gene expression analysis.

    The class represents a "genome" in the form of an ordered set of genes.
    This means that each gene has an index value, i.e. an integer indicating
    its 0-based position in the genome.

    Parameters
    ----------
    genes : Iterable of `ExpGene` objects
        See :attr:`genes` attribute.

    Attributes
    ----------
    genes : list of `ExpGene`
        The genes in the genome. 

    Notes
    -----
    The implementation is very similar to the
    `genometools.basic.GeneSetCollection` class. It uses ordered dictionaries
    to support efficient access by gene name or index, as well as looking up
    the index of specific gene.
    """
    def __init__(self, genes):

        assert isinstance(genes, Iterable)

        self._genes = OrderedDict([g, None] for g in genes)
        for g in self._genes.keys():
            assert isinstance(g, ExpGene)

        """
        # check for duplicate gene names
        counts = Counter(g.name for g in genes)
        for k, v in counts.items():
            if v > 1:
                raise ValueError('Cannot create ExpGenome: more than one gene '
                                 'with name "%s".' % k)
        
        # check for duplicate Ensembl IDs
        counts = Counter(g.ensembl_id for g in self.genes)
        for k, v in counts.items():
            if v > 1 and k is not None:
                raise ValueError('Cannot create ExpGenome: more than one gene '
                                 'with Ensembl ID "%s"' % k)
        """

        # auxiliary variables
        self._genes_by_name = dict([g.name, g] for g in self._genes.keys())
        self._genes_by_index = dict([i, g]
                                    for i, g in
                                    enumerate(self._genes.keys()))
        self._gene_indices = dict([g, i]
                                  for i, g in enumerate(self._genes.keys()))

        logger.debug('Initialized ExpGenome with %d genes.', len(self))

    def __repr__(self):
        return '<%s object (%d genes, hash="%s")>' \
                % (self.__class__.__name__, len(self), self.hash)

    def __str__(self):
        return '<%s object with %d genes>' \
               % (self.__class__.__name__, len(self))

    def __len__(self):
        return len(self._genes)

    def __iter__(self):
        return iter(self._genes.keys())

    def __getitem__(self, key):
        """Simple interface for accessing genes.

        Depending on whether key is an integer or not, look up a gene
        either by index, or by name.
        """
        if isinstance(key, (int, np.integer)):
            if key < 0 or key >= len(self):
                raise ValueError('Index %d out of bounds!' % key)
            return self._genes_by_index[key]
        elif isinstance(key, (str, _oldstr)):
            try:
                return self._genes_by_name[key]
            except KeyError:
                raise ValueError('No gene with name "%s"!' % key)            
        else:
            raise TypeError('Key must be int or str.')

    def __eq__(self, other):
        if self is other:
            return True
        elif type(self) is type(other):
            return repr(self) == repr(other)
        else:
            return NotImplemented

    def __ne__(self, other):
        return not self.__eq__(other)

    def __contains__(self, gene_or_name):
        if isinstance(gene_or_name, (str, _oldstr)):
            return gene_or_name in self._genes_by_name
        elif isinstance(gene_or_name, ExpGene):
            return gene_or_name in self.genes
        else:
            raise TypeError('Must be string or ExpGene instance.')

    @property
    def hash(self):
        """Returns an MD5 hash value for the genome."""
        data_str = ';'.join(repr(g) for g in self.genes)
        data = data_str.encode('UTF-8')
        return str(hashlib.md5(data).hexdigest())

    @property
    def genes(self):
        """Returns a list with all genes."""
        return list(self._genes.keys())

    @property
    def gene_names(self):
        """Returns a list of all gene names."""
        return [g.name for g in self._genes.keys()]

    @property
    def gene_set(self):
        """Returns a set of all genes."""
        return set(self._genes.keys())

    @classmethod
    def from_gene_names(cls, names):
        """Generate a genome from a list of gene names.

        Parameters
        ----------
        names : Iterable of str
            The list of gene names.
        
        Returns
        -------
        ExpGenome
            The genome.
        """
        assert isinstance(names, Iterable)

        genome = cls([ExpGene(n) for n in names])
        return genome


    def index(self, gene_or_name):
        """Returns the index of a given gene.

        The index is 0-based, so the first gene in the genome has the index 0,
        and the last one has index ``len(genome) - 1``.

        Parameters
        ----------
        gene_or_name : str
            The gene or its name (symbol).

        Returns
        -------
        int
            The gene index.
        """
        if isinstance(gene_or_name, (str, _oldstr)):
            # name specified
            try:
                gene = self._genes_by_name[gene_or_name]
            except KeyError:
                raise ValueError('No gene with name "%s"!' % gene_or_name)
        elif isinstance(gene_or_name, ExpGene):
            gene = gene_or_name
        else:
            raise TypeError('Must be a gene name or an ExpGene instance.')

        return self._gene_indices[gene]


    @classmethod
    def read_tsv(cls, path, encoding='UTF-8'):
        """Read genes from tab-delimited text file.

        Parameters
        ----------
        path : str
            The path of the text file.
        encoding : str, optional
            The file encoding. ('UTF-8')

        Returns
        -------
        None
        """
        df = pd.read_csv(path, sep='\t', na_values=['nan', 'NaN'], keep_default_na=False)
        df.columns = [c.lower() for c in df.columns]
        genes = []
        for _, row in df.iterrows():
            genes.append(ExpGene.from_dict(row.to_dict()))

        return cls(genes)


    def write_tsv(self, path, encoding='UTF-8', overwrite=False):
        """Write genes to tab-delimited text file in alphabetical order.

        Parameters
        ----------
        path : str
            The path of the output file.
        encoding: str, optional
            The file encoding. ["UTF-8"]

        Returns
        -------
        None
        """

        # create pandas data frame from the genes 
        data = OrderedDict([i, g.to_dict()]
                           for i, g in enumerate(self._genes.keys()))
        df = pd.DataFrame.from_dict(data, orient='index')

        # write to tab-delimited text file
        sep = '\t'
        if six.PY2:
            sep = sep.encode('UTF-8')

        df.to_csv(path, sep=sep, index=False)

        logger.info('Wrote %d genes to file "%s".', len(self), path)
