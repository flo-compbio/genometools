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

import logging

logger = logging.getLogger(__name__)

class ExpGene(object):
    """A gene in a gene expression analysis.

    Parameters
    ----------
    name: str or unicode
        See :attr:`name` attribute.
    chromosomes: list or tuple of (str or unicode), optional
        See :attr:`chromosomes` attribute.
    ensembl_ids: list or tuple of (str or unicode), optional
        See :attr:`ensembl_ids` attribute.

    Attributes
    ----------
    name: str or unicode
        The gene name (use the official gene symbol, if available).
    chromosomes: None or tuple of (str or unicode)
        The chromosome(s) that the gene is located on.
    ensembl_ids: None or tuple of (str or unicode)
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
    def __init__(self, name, chromosomes = None, ensembl_ids = None):

        if chromosomes is None:
            chromosomes = []
        if ensembl_ids is None:
            ensembl_ids = []

        # checks
        assert isinstance(name, (str, unicode))
        assert isinstance(chromosomes, (list, tuple))
        for chrom in chromosomes:
            assert isinstance(chrom, (str, unicode))

        assert isinstance(ensembl_ids, (list, tuple))
        for id_ in ensembl_ids:
            assert isinstance(id_, (str, unicode))

        self.name = name
        self.chromosomes = tuple(chromosomes)
        self.ensembl_ids = tuple(ensembl_ids)

    def __repr__(self):
        return '<%s "%s" (%s; %s)>' %(self.__class__.__name__,
                self.name, repr(self.chromosomes), repr(self.ensembl_ids))

    def __str__(self):
        chrom_str = 'None'
        if self.chromosomes is not None:
            chrom_str = ','.join(self.chromosomes)
        ens_str = 'None'
        if self.ensembl_ids is not None:
            ens_str = ','.join(self.ensembl_ids)
        return '<%s "%s" (Chromosome(s): %s, EnsemblID(s): %s)>' \
                %(self.__class__.__name__, self.name, chrom_str, ens_str)

    def __hash__(self):
        data = []
        data.append(self.name)
        data.append(self.chromosomes)
        data.append(self.ensembl_ids)

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

    def to_list(self):
        return [self.name, ','.join(self.chromosomes),
                ','.join(self.ensembl_ids)]

    @classmethod
    def from_list(cls, l):
        """Generate an ExpGene object from a list of strings.

        Parameters
        ----------
        l: list or tuple of (str or unicode)
            A list of strings representing gene name, chromosome(s), and
            Ensembl ID(s), respectively. See also :meth:`to_list`.

        Returns
        -------
        `genometools.expression.ExpGene`
        """
        assert isinstance(l, (list, tuple))
        assert len(l) == 3
        for i in range(3):
            assert isinstance(l[i], (str, unicode))

        assert l[0] is not None # name has to be set

        chrom = None
        if l[1] is not None and l[1]:
            chrom = l[1].split(',')

        ens = None
        if l[2] is not None and l[2]:
            ens = l[2].split(',')

        return cls(l[0], chromosomes = chrom, ensembl_ids = ens)
