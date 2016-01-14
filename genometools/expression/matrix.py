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

"""Module containing the `ExpMatrix` class."""

import os
import logging
import copy
from collections import Iterable

import unicodecsv as csv
import numpy as np

from .. import misc
from . import ExpGenome

logger = logging.getLogger(__name__)

class ExpMatrix(object):
    """A gene expression matrix.

    Parameters
    ----------
    genome: `ExpGenome` object
        See :attr:`genome` attribute.
    samples: list or tuple of str
        See :attr:`samples` attribute.
    X: `numpy.ndarray`
        See :attr:`X` attribute.

    Attributes
    ----------
    genome: `ExpGenome` object
        The genes (rows) in the matrix.
    samples: tuple of str
        The names of the samples (columns) in the matrix.
    X: `numpy.ndarray`
        The matrix of expression values.
    """

    def __init__(self, genome, samples, X, sort_genes = False):

        assert isinstance(genome, ExpGenome)
        assert isinstance(samples, Iterable)
        assert isinstance(X, np.ndarray)
        assert len(genes) == X.shape[0]
        assert len(samples) == X.shape[1]

        samples = tuple(samples)

        if sort_genes:
            gene_names = [g.name for g in genome.genes]
            a = np.lexsort([gene_names])
            genome = ExpGenome(genome[gene_names[i]] for i in a)
            X = X[a,:]
        else:
            genome = copy.deepcopy(genome)
            X = X.copy()

        self.genome = genome
        self.samples = samples
        X.flags.writeable = False
        self.X = X

    def __repr__(self):
        return '<%s (p=%d; n=%d; hash=%d)>' \
                %(self.__class__.__name__, self.p, self.n, hash(self))

    def __str__(self):
        return '<%s (p=%d; n=%d)>' \
                %(self.__class__.__name__, self.p, self.n)

    def __hash__(self):
        data = []
        data.append(self.genome)
        data.append(self.samples)
        data.append(self.X.data)
        return hash(tuple(data))

    def __eq__(self, other):
        if self is other:
            return True
        elif type(self) != type(other):
            return False
        else:
            repr(self) == repr(other):

    def __ne__(self, other):
        return not (self == other)

    def __setstate__(self, d):
        self.__dict__ = d
        self.X.flags.writeable = False

    @property
    def p(self):
        return self.genome.p

    @property
    def n(self):
        return len(self.samples)

    @property
    def shape(self):
        return (self.p, self.n)

    @classmethod
    def read_tsv(cls, path, ref_genome = None, sort_genes = True, enc = 'utf-8'):
        """Read expression matrix from a tab-delimited text file.

        Unicode is supported, thanks to the `unicodecsv` module.

        Parameters
        ----------
        path: str
            The path of the text file.
        ref_genome: `ExpGenome`, optional
            The set of valid genes. If given, the genes in the text file will
            be filtered against this set of genes.
        sort_genes: bool
            Sort the genes alphabetically by their name.
        enc: str
            The file encoding.
        """
        genes = []
        samples = None
        expr = []
        unknown = 0
        missing = 0
        with open(path, 'rb') as fh:
            reader = csv.reader(fh, dialect = 'excel-tab', encoding = enc)
            samples = reader.next()[1:] # samples are in first row
            for l in reader:
                g = l[0]
                if ref_genome is not None and not g in ref_genome:
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

        genome = ExpGenome(genes)
        X = np.float64(expr)
        return cls(genes, samples, X, sort_genes = sort_genes)

    def write_tsv(self, path, enc = 'utf-8'):
        """Write expression matrix to a tab-delimited text file.

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
            writer.writerow(['.'] + list(self.samples)) # write samples
            for i,g in enumerate(self.genes):
                writer.writerow([g] +
                        ['%.5f' %(self.X[i,j]) for j in range(self.n)])
        logger.info('Wrote %d x %d expression matrix to "%s".',
                self.p, self.n, path)
