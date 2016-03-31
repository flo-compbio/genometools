# Copyright (c) 2016 Florian Wagner
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

"""Module containing the `ExpProfile` class."""

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
matrix = importlib.import_module('.matrix', package='genometools.expression')
# "from . import matrix" does not work, due to cyclical imports

logger = logging.getLogger(__name__)

class ExpProfile(pd.Series):
    """A gene expression profile.

    Parameters
    ----------
    x: 1-dimensional `numpy.ndarray`
        See :attr:`x` attribute.
        
    Keyword-only Parameters
    -----------------------
    genes: list or tuple of str
        See :attr:`genes` attribute.
    name: str
        See :attr:`name` attribute.
        
    Additional Parameters
    ---------------------
    All `pandas.Series` parameters.

    Attributes
    ----------
    x: 1-dimensional `numpy.ndarray`
        The vector with expression values.
    genes: tuple of str
        The names of the genes (rows) in the matrix.
    label: str
        The sample label.
    """
    def __init__(self, *args, **kwargs):
        
        # check if user provided "x" keyword argument
        x = kwargs.pop('x', None)
        if x is not None:
            assert isinstance(x, np.ndarray)
            kwargs['data'] = x
        
        # check if user provided "genes" keyword argument
        genes = kwargs.pop('genes', None)
        if genes is not None:
            assert isinstance(genes, (list, tuple))
            for g in genes:
                assert isinstance(g, str)
        
        # check if user provided "label" keyword argument
        label = kwargs.pop('label', None)
        if label is not None:
            assert isinstance(label, str)
         
        # call base class constructor
        pd.Series.__init__(self, *args, **kwargs)
        
        if genes is not None:
            # set (overwrite) index with user-provided list
            self.index = genes

        if label is not None:
            # set (overwrite) series name with user-provided sample label
            self.name = label
                        
    def __hash__(self):
        # warning: involves copying all the data
        data = []
        data.append(tuple(self.genes))
        data.append(tuple(self.label))
        x = self.x.copy()
        x.flags.writeable = False
        data.append(x.data)
        return hash(tuple(data))

    @property
    def _constructor(self):
        return ExpProfile

    @property
    def _constructor_expanddim(self):
        return matrix.ExpMatrix
    
    @property
    def p(self):
        """The number of genes."""
        #return len(self.genes)
        return self.shape[0]

    @property
    def genes(self):
        """Returns the gene names (= series index) as a list."""
        return self.index.values.tolist()

    @genes.setter
    def genes(self, genes):
        self.index = genes

    @property
    def label(self):
        """Returns the sample label (= series name)."""
        return self.name

    @label.setter
    def label(self, label):
        self.name = label

    @property
    def x(self):
        """Returns the expression values as a numpy vector."""
        return self.values

    @x.setter
    def x(self, x):
        self.x[:] = x

    def get_genome(self):
        """Get an ExpGenome representation of the genes in the profile.

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
        """Sort the rows of the profile alphabetically by gene name.

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
        """Read expression profile from a tab-delimited text file.

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
        `ExpProfile`
            The expression profile.
        """
        # checks
        assert isinstance(path, str)
        if genome is not None:
            assert isinstance(genome, ExpGenome)
        assert isinstance(encoding, str)

        # "squeeze = True" ensures that a pd.read_tsv returns a series
        # as long as there is only one column
        e = cls(pd.read_csv(path, sep = '\t', index_col = 0, header = 0,
                         encoding = encoding, squeeze = True))

        if genome is not None:
            # filter genes
            e = e.filter_against_genome(genome)

        return e

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

        sep = '\t'
        if six.PY2:
            sep = sep.encode('UTF-8')

        self.to_csv(
            path, sep = sep, float_format = '%.5f', mode = 'w',
            encoding = encoding, header = True
        )

        logger.info('Wrote expression profile "%s" with %d genes to "%s".',
                self.name, self.p, path)
