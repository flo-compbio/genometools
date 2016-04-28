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
from builtins import *

import os
import logging
import hashlib
from collections import OrderedDict

import unicodecsv as csv
import numpy as np

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
    def __init__(self, exp_genes):

        assert isinstance(exp_genes, (list, tuple))
        for g in exp_genes:
            assert isinstance(g, ExpGene)

        all_genes = [eg.name for eg in exp_genes]
        if len(all_genes) != len(set(all_genes)):
            raise ValueError('Cannot create ExpGenome: Not all genes have a '
                             'unique name.')

        self._exp_gene_dict = OrderedDict([eg.name, eg] for eg in exp_genes)
        self._gene_names = tuple(self._exp_gene_dict.keys())
        self._gene_indices = OrderedDict([eg.name, i]
                                         for i, eg in enumerate(exp_genes))
        logger.debug('Initialized ExpGenome with %d genes.', len(self))

    def __repr__(self):
        return '<%s object (%d genes; hash=%s)>' \
                % (self.__class__.__name__, len(self), self.hash)

    def __str__(self):
        return '<%s object with %d genes>' \
               % (self.__class__.__name__, len(self))

    def __len__(self):
        return len(self._exp_gene_dict)

    def __iter__(self):
        return iter(self._exp_gene_dict.values())

    def __getitem__(self, key):
        """Simple interface for accessing genes.

        Depending on whether key is an integer or not, look up a gene
        either by index, or by name.
        """
        if isinstance(key, (int, np.integer)):
            return self.get_by_index(key)
        else:
            return self.get_by_name(key)

    @property
    def hash(self):
        data_str = ';'.join(repr(g) for g in self._exp_gene_dict)
        data = data_str.encode('ascii')
        return hashlib.md5(data).hexdigest()

    def __eq__(self, other):
        if self is other:
            return True
        elif type(self) is type(other):
            return self.__dict__ == other.__dict__
        else:
            return NotImplemented

    def __ne__(self, other):
        return not (self == other)

    def __contains__(self, gene):
        return gene.name in self._exp_gene_dict

    @property
    def exp_genes(self):
        """Returns a list with all genes."""
        return list(self._exp_gene_dict.values())

    @property
    def genes(self):
        """Returns a list with all gene names."""
        return [eg.name for eg in self._exp_gene_dict.values()]

    @classmethod
    def from_gene_names(cls, genes):
        """Generate a genome from a list of gene names.

        Parameters
        ----------
        genes: list or tuple of str
            The list of gene names.
        
        Returns
        -------
        ExpGenome
            The genome.
        """
        assert isinstance(genes, (list, tuple))
        for g in genes:
            assert isinstance(g, str)

        genome = cls([ExpGene(g) for g in genes])
        return genome

    def get_by_name(self, gene):
        """Look up a gene by its name.

        Parameters
        ----------
        gene: str
            The gene name.
        
        Returns
        -------
        `genometools.expression.ExpGene`
            The gene.
        """
        if not isinstance(gene, str):
            raise ValueError('Gene name must be a string.')

        try:
            return self._exp_gene_dict[gene]
        except KeyError:
            raise ValueError('No gene with name "%s"!' % gene)

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
        if not isinstance(i, (int, np.integer)):
            raise ValueError('Index must be an integer.')

        if i >= len(self):
            raise ValueError('Index %d out of bounds for genome with %d genes.'
                             % (i, len(self)))

        return self._exp_gene_dict[self._gene_names[i]]

    def index(self, gene):
        """Returns the index of the gene with the given gene.

        The index is 0-based, so the first gene in the genome has the index 0,
        and the last one has index p-1.

        Parameters
        ----------
        gene: str
            The gene.

        Returns
        -------
        int
            The gene index.
        """
        if not isinstance(gene, str):
            raise ValueError('Gene name must be a string.')

        try:
            return self._gene_indices[gene]
        except KeyError:
            raise ValueError('No gene with name "%s"!' % gene)

    @classmethod
    def read_tsv(cls, path, encoding='UTF-8'):
        """Read genes from tab-delimited text file.

        Parameters
        ----------
        path: str
            The path of the text file.
        encoding: str, optional
            The file encoding. ("UTF-8")

        Returns
        -------
        None
        """
        genes = []
        with open(path, 'rb') as fh:
            reader = csv.reader(fh, dialect='excel-tab', encoding=encoding)
            for l in reader:
                genes.append(ExpGene.from_list(l))
        return cls(genes)

    def write_tsv(self, output_file, encoding='UTF-8'):
        """Write genes to tab-delimited text file in alphabetical order.

        Parameters
        ----------
        output_file: str
            The path of the output file.
        encoding: str, optional
            The file encoding. ("UTF-8")

        Returns
        -------
        None
        """
        with open(output_file, 'wb') as ofh:
            writer = csv.writer(
                ofh, dialect='excel-tab', encoding=encoding,
                lineterminator=os.linesep, quoting=csv.QUOTE_NONE
            )
            for eg in self._exp_gene_dict.values():
                writer.writerow(eg.to_list())
        logger.info('Wrote %d genes to file "%s".', len(self), output_file)
