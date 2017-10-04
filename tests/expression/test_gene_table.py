# Copyright (c) 2016, 2017 Florian Wagner
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

"""Tests for the `ExpGeneTable` class."""

from copy import deepcopy
from collections import OrderedDict

import pytest

from genometools.misc import get_logger
from genometools.expression import ExpGene, ExpGeneTable

# _LOGGER = get_logger(__name__, verbose=True)


def test_init(my_gene_table, my_genes):
    assert isinstance(my_gene_table, ExpGeneTable)
    assert isinstance(repr(my_gene_table), str)
    assert isinstance(str(my_gene_table), str)
    assert isinstance(my_gene_table.hash, str)
    assert len(my_gene_table) == len(my_genes)

    assert my_gene_table.genes == my_genes

    other = deepcopy(my_gene_table)
    assert other is not my_gene_table
    assert other.equals(my_gene_table)


def test_from_ids_and_names(my_gene_ids, my_gene_names):
    gene_names = OrderedDict([id_, name]
                      for id_, name in zip(my_gene_ids, my_gene_names))
    gene_table = ExpGeneTable.from_gene_ids_and_names(gene_names)
    assert len(gene_table) == len(my_gene_names)
    for id_, name, g in zip(my_gene_ids, my_gene_names, gene_table.genes):
        assert g.ensembl_id == id_
        assert g.name == name
        assert g.chromosome is None
        assert g.position is None
        assert g.length is None

def test_access(my_gene_table, my_genes, my_unknown_gene_name):
    for i, g in enumerate(my_genes):
        assert g in my_gene_table
        assert my_gene_table.genes[i] == g
        #assert my_genome[g.name] is g
        #assert my_genome.index(g) == i

    #with pytest.raises(ValueError):
    #    # test invalid gene name
    #    my_gene_table.genes[my_unknown_gene_name]

    #with pytest.raises(ValueError):
    #    # test invalid gene index
    #    my_gene_table[len(my_genes)]

    #with pytest.raises(ValueError):
    #    # test getting index of a gene that does not exist
    #    my_gene_table.index(my_unknown_gene_name)


def test_tsv(tmpdir, my_gene_table):
    tmp_file = str(tmpdir.join('_genome.tsv'))
    # print(type(_genome.exp_genes[0]))
    my_gene_table.write_tsv(tmp_file)
    other = ExpGeneTable.read_tsv(tmp_file)
    assert my_gene_table.equals(other)