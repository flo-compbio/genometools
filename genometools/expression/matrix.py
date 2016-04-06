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

"""Module containing the `ExpMatrix` class."""

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import *

import os
import io
import logging
import copy
import importlib

import pandas as pd
import numpy as np
import unicodecsv as csv
import six

from .. import misc
from . import ExpGene, ExpGenome
profile = importlib.import_module('.profile', package='genometools.expression')
# - "import profile" is not possible, since a "profile" module exists
#    in the standard library
# - "from . import profile" fails due to cyclical imports

logger = logging.getLogger(__name__)

class ExpMatrix(pd.DataFrame):
    """A gene expression matrix.

    Parameters
    ----------
    X: 2-dimensional `numpy.ndarray`
        See :attr:`X` attribute.

    Keyword-only Parameters
    -----------------------
    genes: list or tuple of str
        See :attr:`genes` attribute.
    samples: list or tuple of str
        See :attr:`samples` attribute.

    Additional Parameters
    ---------------------
    All `pandas.DataFrame` parameters.

    Attributes
    ----------
    genes: tuple of str
        The names of the genes (rows) in the matrix.
    samples: tuple of str
        The names of the samples (columns) in the matrix.
    X: 2-dimensional `numpy.ndarray`
        The matrix of expression values.
    """
    def __init__(self, *args, **kwargs):
        
        # check if user provided "X" keyword argument
        X = kwargs.pop('X', None)
        if X is not None:
            assert isinstance(X, np.ndarray)
            kwargs['data'] = X

        # check if user provided "genes" keyword argument
        genes = kwargs.pop('genes', None)
        if genes is not None:
            assert isinstance(genes, (list, tuple))
            for g in genes:
                assert isinstance(g, str)
            
        # check if user provided "samples" keyword argument
        samples = kwargs.pop('samples', None)
        if samples is not None:
            assert isinstance(samples, (list, tuple))
            for s in samples:
                assert isinstance(s, str)

        
        if genes is not None:
            kwargs['index'] = genes

        if samples is not None:
            kwargs['columns'] = samples

        # call base class constructor
        pd.DataFrame.__init__(self, *args, **kwargs)
        
        #if genes is not None:
        #    # set (overwrite) index with user-provided list
        #    self.index = genes
            
        #if samples is not None:
        #    # set (overwrite) index with user-provided list
        #    self.columns = samples
        
    def __hash__(self):
        # warning: involves copying all the data
        data = []
        data.append(tuple(self.genes))
        data.append(tuple(self.samples))
        data.append(self.X.tobytes())
        return hash(tuple(data))

    @property
    def _constructor(self):
        return ExpMatrix
    
    @property
    def _constructor_sliced(self):
        return profile.ExpProfile

    @property
    def p(self):
        """The number of genes."""
        return self.shape[0]

    @property
    def n(self):
        """The number of samples."""
        return self.shape[1]

    @property
    def genes(self):
        """Returns the gene (row) names as a list."""
        return self.index.values.tolist()

    @genes.setter
    def genes(self, gene_list):
        self.index = gene_list

    @property
    def samples(self):
        """Returns the sample (column) names as a list."""
        return self.columns.values.tolist()

    @samples.setter
    def samples(self, sample_list):
        self.columns = sample_list

    @property
    def X(self):
        """Returns the expression values as a numpy array."""
        return self.values

    @X.setter
    def X(self, X):
        self.values[:,:] = X

    def get_genome(self):
        """Get a ExpGenome representation of the genes in the matrix.

        Parameters
        ----------
        None

        Returns
        -------
        `genometools.expression.ExpGenome`
            The genome.
        """

        genes = [ExpGene(g) for g in self.genes]
        genome = ExpGenome(genes)
        return genome

    def sort_genes(self, stable = False):
        """Sort the rows of the matrix alphabetically by gene name.

        Parameters
        ----------
        stable: bool, optional
            If set to True, uses a stable sorting algorithm. (False)
        
        Returns
        -------
        None
        """
        kind = 'quicksort'
        if stable:
            kind = 'mergesort'
        self.sort_index(kind = kind)

    def center_genes(self):
        """Center the expression of each gene (row)."""
        self.X = self.X - np.tile(np.mean(self.X, axis = 1), (self.n, 1)).T

    def standardize_genes(self):
        """Standardize the expression of each gene (row)."""
        self.center_genes()
        self.X = self.X / np.tile(np.std(self.X, axis = 1, ddof = 1), 
                                  (self.n, 1)).T

    def filter_against_genome(self, genome):
        """Filter the expression matrix against a genome (set of genes).

        Parameters
        ----------
        genome: `genometools.expression.ExpGenome`
            The genome to filter the genes against.

        Returns
        -------
        ExpMatrix
            The filtered expression matrix.
        """
        assert isinstance(genome, ExpGenome)

        return self.loc[self.index & genome.genes]

    @classmethod
    def read_tsv(cls, path, genome = None, encoding = 'UTF-8'):
        """Read expression matrix from a tab-delimited text file.

        Parameters
        ----------
        path: str
            The path of the text file.
        genome: `ExpGenome` object, optional
            The set of valid genes. If given, the genes in the text file will
            be filtered against this set of genes. (None)
        encoding: str, optional
            The file encoding. ("UTF-8")

        Returns
        -------
        `ExpMatrix`
            The expression matrix.
        """
        # checks
        assert isinstance(path, str)
        if genome is not None:
            assert isinstance(genome, ExpGenome)
        assert isinstance(encoding, str)

        # use pd.read_csv to parse the tsv file into a DataFrame
        E = cls(pd.read_csv(path, sep = '\t', index_col = 0, header = 0, encoding = encoding))

        if genome is not None:
            # filter genes
            E = E.filter_against_genome(genome)

        return E

    def write_tsv(self, path, encoding = 'UTF-8'):
        """Write expression matrix to a tab-delimited text file.

        Parameters
        ----------
        path: str
            The path of the output file.
        encoding: str, optional
            The file encoding. ("UTF-8")

        Returns
        -------
        None
        """
        assert isinstance(path, str)
        assert isinstance(encoding, str)

        #sep = str('\t')
        sep = '\t'
        if six.PY2:
            sep = sep.encode('UTF-8')

        self.to_csv(
            path, sep = sep, float_format = '%.5f', mode = 'w',
            encoding = encoding, quoting = csv.QUOTE_NONE
        )

        logger.info('Wrote %d x %d expression matrix to "%s".',
                self.p, self.n, path)
