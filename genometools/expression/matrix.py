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

import csv
import logging

import numpy as np

logger = logging.getLogger(__name__)

class ExpMatrix(object):
    """A gene expression matrix.

    Parameters
    ----------
    genes: list or tuple of str
        See :attr:`genes` attribute.
    samples: list or tuple of str
        See :attr:`samples` attribute.
    E: `numpy.ndarray`
        See :attr:`E` attribute.

    Attributes
    ----------
    genes: tuple of str
        The names of the genes (rows) in the matrix.
    samples: tuple of str
        The names of the samples (columns) in the matrix.
    E: `numpy.ndarray`
        The matrix of expression values.
    """

    def __init__(self,genes,samples,E,preserve_gene_order=None):

        assert len(genes) == E.shape[0]
        assert len(samples) == E.shape[1]

        genes = tuple(genes)
        samples = tuple(samples)
        E = E.copy()

        if preserve_gene_order is None or (not preserve_gene_order):
            # make sure genes are in alphabetical order
            a = np.lexsort([genes])
            if not np.all(a == np.arange(len(genes))):
                genes = tuple(genes[i] for i in a)
                E = E[a,:]

        self.genes = genes
        self.samples = samples
        self.E = E

    @property
    def p(self):
        return len(self.genes)

    @property
    def n(self):
        return len(self.samples)

    @property
    def shape(self):
        return (len(self.genes),len(self.samples))

    def __repr__(self):
        self.E.flags.writable = False
        h = hash((self.genes,self.samples,self.E))
        self.E.flags.writable = True
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
    def read_tsv(cls, path, genome = None, preserve_gene_order = None):
        """Read expression matrix from a tab-delimited text file.

        Parameters
        ----------
        path: str
            The path of the text file.
        genome: `ExpGenome`, optional
            The set of valid genes. If given, the genes in the text file will
            be filtered against the set of genes in the genome.
        """
        genes = []
        samples = None
        expr = []
        unknown = 0
        missing = 0
        with open(path) as fh:
            reader = csv.reader(fh,dialect='excel-tab')
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

        E = np.float64(expr)
        return cls(genes, samples, E,
                preserve_gene_order = preserve_gene_order)

    def write_tsv(self, path):
        """Write expression matrix to a tab-delimited text file.

        Parameters
        ----------
        path: str
            The path of the output file.
        """
        with open(path,'w') as ofh:
            writer = csv.writer(ofh, dialect = 'excel-tab',
                    lineterminator = '\n', quoting = csv.QUOTE_NONE)
            writer.writerow(['.'] + list(self.samples)) # write samples
            for i,g in enumerate(self.genes):
                writer.writerow([g] +
                        ['%.5f' %(self.E[i,j]) for j in range(self.n)])
        logger.info('Wrote %d x %d expression matrix to "%s".',
                self.p, self.n, path)
