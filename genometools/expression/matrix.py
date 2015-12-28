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

import unicodecsv as csv
import numpy as np

from genometools import misc

logger = logging.getLogger(__name__)

class ExpMatrix(object):
    """A gene expression matrix.

    Parameters
    ----------
    genes: list or tuple of str
        See :attr:`genes` attribute.
    samples: list or tuple of str
        See :attr:`samples` attribute.
    X: `numpy.ndarray`
        See :attr:`X` attribute.

    Attributes
    ----------
    genes: tuple of str
        The names of the genes (rows) in the matrix.
    samples: tuple of str
        The names of the samples (columns) in the matrix.
    X: `numpy.ndarray`
        The matrix of expression values.
    """

    def __init__(self, genes, samples, X, preserve_gene_order = False):

        assert len(genes) == X.shape[0]
        assert len(samples) == X.shape[1]

        genes = tuple(genes)
        samples = tuple(samples)
        X = X.copy()

        if (not preserve_gene_order):
            # make sure genes are in alphabetical order
            a = np.lexsort([genes])
            if not np.all(a == np.arange(len(genes))):
                genes = tuple(genes[i] for i in a)
                X = X[a,:]

        self.genes = genes
        self.samples = samples
        self.X = X

    @property
    def p(self):
        return len(self.genes)

    @property
    def n(self):
        return len(self.samples)

    @property
    def shape(self):
        return (len(self.genes), len(self.samples))

    def __repr__(self):
        self.X.flags.writeable = False
        h = hash((self.genes, self.samples, self.X.data))
        self.X.flags.writeable = True
        return '<ExprMatrix with %d genes, %d samples (hash = %d)>' \
                %(self.p,self.n,h)

    def __str__(self):
        return '<ExprMatrix with %d genes, %d samples>' \
                %(self.p,self.n)

    def __hash__(self):
        return hash(repr(self))

    def __eq__(self,other):
        if type(self) != type(other):
            return False
        elif repr(self) == repr(other):
            return True
        else:
            return False

    @classmethod
    def read_tsv(cls, path, genome = None, preserve_gene_order = False, enc = 'utf-8'):
        """Read expression matrix from a tab-delimited text file.

        Unicode is supported, thanks to the `unicodecsv` module.

        Parameters
        ----------
        path: str
            The path of the text file.
        genome: `ExpGenome`, optional
            The set of valid genes. If given, the genes in the text file will
            be filtered against the set of genes in the genome.
        preserve_gene_order: bool
            Do not sort the genes alphabetically.
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
                if genome is not None and not genome.has_gene(g):
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
        return cls(genes, samples, X,
                preserve_gene_order = preserve_gene_order)

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
