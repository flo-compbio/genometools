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

"""Module containing the `ExpGene` class.

"""

import logging

logger = logging.getLogger(__name__)

class ExpGene(object):
    """A gene in a gene expression analysis.

    Parameters
    ----------
    name: str
        See :attr:`name` attribute.
    chromosomes: list of str, optional
        See :attr:`chromosomes` attribute.
    ensembl_ids: list of str, optional
        See :attr:`ensembl_ids` attribute.

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
    The reason :attr:`chromosomes` and :attr:`ensembl_ids` are lists is mainly
    to accommodate genes located on the pseudoautosomal region of the X/Y
    chromosomes. Although these genes have separate Ensembl IDs depending
    for their X and Y "versions", in gene expression analyses they should be
    treated as the same gene. This class therefore represents a more "abstract"
    idea of a gene, not its physical manifestation in the genome.
    """
    def __init__(self, name, chromosomes = None, ensembl_ids = None):

        self.name = name

        if chromosomes is None:
            chromosomes = []
        self.chromosomes = chromosomes
        if ensembl_ids is None:
            ensembl_ids = []
        self.ensembl_ids = ensembl_ids

    def __repr__(self):
        return '<ExpGene "%s" (%s; %s)>' \
                %(self.name, repr(self.chromosomes), repr(self.ensembl_ids))

    def __str__(self):
        chrom_str = 'None'
        if self.chromosomes is not None:
            chrom_str = ','.join(self.chromosomes)
        ens_str = 'None'
        self.ensembl_ids is not None:
            ens_str = ','.join(self.ensembl_ids)
        return '<ExpGene "%s" (Chromosome(s): %s, EnsemblID(s): %s)>' \
                %(self.name, chrom_str, ens_str)

    def __hash__(self):
        data = []
        data.append(self.name)

        if self.chromosomes is None:
            data.append(None)
        else:
            data.append(tuple(self.chromosomes))

        if self.ensembl_ids is None:
            data.append(None)
        else:
            data.append(tuple(self.ensembl_ids))

        return hash(tuple(data))

    def __eq__(self, other):
        if self is other:
            return True
        elif type(self) != type(other):
            return False
        else:
            return repr(self) == repr(other)

    def to_list(self):
        return [self.name, ','.join(self.chromosomes),
                ','.join(self.ensembl_ids)]

    @classmethod
    def from_list(cls, l):
        assert len(l) == 3
        assert l[0]

        chrom = None
        if l[1]:
            chrom = l[1].split(',')

        ens = None
        if l[2]:
            ens = l[2].split(',')

        return cls(l[0], chromosomes = chrom, ensembl_ids = ens)
