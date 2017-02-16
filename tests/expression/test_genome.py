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

"""Tests for the `ExpGenome` class."""

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import str as text

from copy import deepcopy
from collections import Iterable

import pytest

from genometools.misc import get_logger
from genometools.expression import ExpGene, ExpGenome

logger = get_logger(__name__, verbose=True)


def test_init(my_genome, my_genes):
    assert isinstance(my_genome, ExpGenome)
    assert isinstance(repr(my_genome), str)
    assert isinstance(str(my_genome), str)
    assert isinstance(text(my_genome), text)
    assert isinstance(my_genome.hash, text)
    assert len(my_genome) == len(my_genes)
    assert isinstance(my_genome, Iterable)

    assert my_genome.genes == my_genes

    other = deepcopy(my_genome)
    assert other is not my_genome
    assert other == my_genome


def test_from_names(my_gene_names):
    genome = ExpGenome.from_gene_names(my_gene_names)
    assert len(genome) == len(my_gene_names)
    for i, (n, g) in enumerate(zip(my_gene_names, genome)):
        assert g.name == n
        assert g.chromosome is None
        assert g.position is None
        assert g.length is None
        assert g.ensembl_id is None

def test_access(my_genome, my_genes, my_unknown_gene_name):
    for i, g in enumerate(my_genes):
        assert my_genome[i] is g
        assert my_genome[g.name] is g
        assert my_genome.index(g) == i
        assert g in my_genome

    for i, g in enumerate(my_genome):
        assert g == my_genes[i]

    with pytest.raises(ValueError):
        # test invalid gene name
        my_genome[my_unknown_gene_name]

    with pytest.raises(ValueError):
        # test invalid gene index
        my_genome[len(my_genes)]

    with pytest.raises(ValueError):
        # test getting index of a gene that does not exist
        my_genome.index(my_unknown_gene_name)


def test_tsv(tmpdir, my_genome):
    tmp_file = str(tmpdir.join('_genome.tsv'))
    # print(type(_genome.exp_genes[0]))
    my_genome.write_tsv(tmp_file)
    other = ExpGenome.read_tsv(tmp_file)
    assert my_genome == other