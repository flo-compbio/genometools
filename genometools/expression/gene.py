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

"""Module containing the `ExpGene` class.

"""

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
_oldstr = str
from builtins import *

import logging

logger = logging.getLogger(__name__)


class ExpGene(object):
    """A gene in a gene expression analysis.

    Instances are to be treated as immutable, to allow use of ExpGene
    objects to be used in sets etc.

    Parameters
    ----------
    name: str
        See :attr:`name` attribute.
    chromosomes: list or tuple of str, optional
        See :attr:`chromosomes` attribute. [empty list]
    ensembl_ids: list or tuple of str, optional
        See :attr:`ensembl_ids` attribute. [empty list]

    Attributes
    ----------
    name: str
        The gene name (use the official gene symbol, if available).
    chromosomes: list of str
        The chromosome(s) that the gene is located on.
    ensembl_ids: list of str
        The Ensembl ID(s) of the gene.

    Notes
    -----
    The reason :attr:`chromosomes` and :attr:`ensembl_ids` are lists / tuples
    is mainly to accommodate genes located on the pseudoautosomal region of
    the X/Y chromosomes. Although these genes have separate Ensembl IDs
    for their X and Y "versions", in gene expression analyses they should be
    treated as the same gene. This class therefore represents a more "abstract"
    idea of a gene, not its physical manifestation in the genome.
    """
    def __init__(self, name, chromosomes=None, ensembl_ids=None):

        if chromosomes is None:
            chromosomes = []
        if ensembl_ids is None:
            ensembl_ids = []

        # type checks
        assert isinstance(name, (str, _oldstr))
        assert isinstance(chromosomes, (list, tuple))
        for chrom in chromosomes:
            assert isinstance(chrom, (str, _oldstr))
        assert isinstance(ensembl_ids, (list, tuple))
        for id_ in ensembl_ids:
            assert isinstance(id_, (str, _oldstr))

        self._name = name
        self._chromosomes = tuple(chromosomes)
        self._ensembl_ids = tuple(ensembl_ids)

    def __repr__(self):
        return '<%s instance (name="%s", chromosomes=%s, ensembl_ids=%s)>' \
               % (self.__class__.__name__, self._name,
                  repr(self._chromosomes), repr(self._ensembl_ids))

    def __str__(self):
        return '<%s instance "%s" (Chromosome(s): %s, EnsemblID(s): %s)>' \
               % (self.__class__.__name__, self._name,
                  str(self._chromosomes), str(self._ensembl_ids))

    def __eq__(self, other):
        if self is other:
            return True
        elif type(self) is type(other):
            return repr(self) == repr(other)
        else:
            return NotImplemented

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return hash(repr(self))

    @property
    def name(self):
        return self._name

    @property
    def chromosomes(self):
        return list(self._chromosomes)

    @property
    def ensembl_ids(self):
        return list(self._ensembl_ids)

    def to_list(self):
        return [self._name, ','.join(self._chromosomes),
                ','.join(self._ensembl_ids)]

    @classmethod
    def from_list(cls, data):
        """Generate an ExpGene object from a list of strings.

        Parameters
        ----------
        data: list or tuple of str
            A list of strings representing gene name, chromosome(s), and
            Ensembl ID(s), respectively. See also :meth:`to_list`.

        Returns
        -------
        `ExpGene`
        """
        assert isinstance(data, (list, tuple))
        assert len(data) == 3
        for l in data:
            assert isinstance(l, str)

        chrom = None
        if data[1]:
            chrom = data[1].split(',')

        ens = None
        if data[2]:
            ens = data[2].split(',')

        return cls(data[0], chromosomes=chrom, ensembl_ids=ens)
