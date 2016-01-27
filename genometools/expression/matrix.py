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

import os
import io
import logging
import copy

import unicodecsv as csv
import numpy as np

from .. import misc

logger = logging.getLogger(__name__)

class ExpMatrix(object):
    """A gene expression matrix.

    Parameters
    ----------
    genes: list or tuple of (str or unicode)
        See :attr:`genes` attribute.
    samples: list or tuple of (str or unicode)
        See :attr:`samples` attribute.
    X: 2-dimensional `numpy.ndarray`
        See :attr:`X` attribute.

    Attributes
    ----------
    genes: tuple of (str or unicode)
        The names of the genes (rows) in the matrix.
    samples: tuple of (str or unicode)
        The names of the samples (columns) in the matrix.
    X: 2-dimensional `numpy.ndarray`
        The matrix of expression values.
    """
    def __init__(self, genes, samples, X):

        assert isinstance(genes, (list, tuple))
        for g in genes:
            assert isinstance(g, (str, unicode))
        assert isinstance(samples, (list, tuple))
        for s in samples:
            assert isinstance(s, (str, unicode))
        assert isinstance(X, np.ndarray)
        assert len(genes) == X.shape[0]
        assert len(samples) == X.shape[1]

        self.genes = tuple(genes)
        self.samples = tuple(samples)
        self.X = X.copy()
        self.X.flags.writeable = False

    def __repr__(self):
        return '<%s (p=%d; n=%d; hash=%d)>' \
                %(self.__class__.__name__, self.p, self.n, hash(self))

    def __str__(self):
        return '<%s (p=%d; n=%d)>' \
                %(self.__class__.__name__, self.p, self.n)

    def __hash__(self):
        data = []
        data.append(self.genes)
        data.append(self.samples)
        data.append(self.X.data)
        return hash(tuple(data))

    def __eq__(self, other):
        if self is other:
            return True
        elif type(self) != type(other):
            return False
        else:
            return repr(self) == repr(other)

    def __ne__(self, other):
        return not (self == other)

    def __setstate__(self, d):
        self.__dict__ = d
        # values of ndarray flags are not stored in pickle
        self.X.flags.writeable = False

    @property
    def p(self):
        """The number of genes."""
        return len(self.genes)

    @property
    def n(self):
        """The number of samples."""
        return len(self.samples)

    @property
    def shape(self):
        """The shape of the matrix (# genes, # samples)."""
        return (self.p, self.n)

    def sort(self):
        """Sort the rows of the matrix alphabeticaly by gene name.

        Parameters
        ----------
        None
        
        Returns
        -------
        None
        """
        a = np.lexsort([self.genes])
        self.genes = tuple(self.genes[i] for i in a)
        self.X = self.X[a,:]

    @classmethod
    def read_tsv(cls, path, genome = None, sort_genes = False, enc = 'UTF-8'):
        """Read expression matrix from a tab-delimited text file.

        Unicode is supported, thanks to the `unicodecsv` module.

        Parameters
        ----------
        path: str
            The path of the text file.
        genome: `ExpGenome` object, optional
            The set of valid genes. If given, the genes in the text file will
            be filtered against this set of genes. (None)
        sort_genes: bool, optional
            Also sort the genes alphabetically by their name. (False)
        enc: str
            The file encoding. ("UTF-8")

        Returns
        -------
        `genometools.expression.ExpMatrix`
            The expression matrix.
        """
        genes = []
        samples = None
        expr = []
        unknown = 0
        missing = 0
        with io.open(path, 'rb') as fh:
            reader = csv.reader(fh, dialect = 'excel-tab', encoding = enc)
            samples = reader.next()[1:] # samples are in first row
            for l in reader:
                g = l[0]
                if genome is not None and not g in genome:
                    unknown += 1
                    continue
                if 'NA' in l[1:]:
                    missing += 1
                    continue
                genes.append(g)
                expr.append(l[1:])

        if unknown > 0:
            logger.warning('Ignored %d / %d genes with unknown name (%.1f%%).',
                    unknown, len(genes), 100*(unknown/float(len(genes))))

        if missing > 0:
            logger.warning('Ignored %d / %d genes with missing data (%.1f%%).',
                    missing, len(genes), 100*(missing/float(len(genes))))

        X = np.float64(expr)
        logger.debug('Expression matrix shape: %s', str(X.shape))
        return cls(genes, samples, X)

    def write_tsv(self, path, enc = 'UTF-8'):
        """Write expression matrix to a tab-delimited text file.

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
            writer.writerow(['.'] + list(self.samples)) # write samples
            for i,g in enumerate(self.genes):
                writer.writerow([g] +
                        ['%.5f' %(self.X[i,j]) for j in range(self.n)])
        logger.info('Wrote %d x %d expression matrix to "%s".',
                self.p, self.n, path)
