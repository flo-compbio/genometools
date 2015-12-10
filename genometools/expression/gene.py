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
    """A gene in gene expression analysis.

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
        return '<ExpGene "%s">' %(self.name)

    def __str__(self):
        chrom_str = ','.join(self.chromosomes)
        ensembl_str = ','.join(self.ensembl_ids)
        return '<ExpGene "%s" (Chromosome(s): %s, EnsemblID(s): %s)>' \
                %(self.name, chrom_str, ensembl_str)

    def __hash__(self):
        return hash(repr(self))

    def __eq__(self,other):
        if type(self) != type(other):
            return False
        elif repr(self) == repr(other):
            return True
        else:
            return False

    def to_list(self):
        return [self.name,','.join(self.chromosomes),','.join(self.ensembl_ids)]

    @classmethod
    def from_list(cls,l):
        return cls(l[0],l[1].split(','),l[2].split(','))
