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

"""Tests for the `GeneSetCollection` class."""

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import str as text

from copy import deepcopy

import pytest
from genometools.basic import GeneSet, GeneSetCollection


@pytest.fixture
def my_gene_sets(my_gene_set):
    gene_sets = []
    for i in range(3):
        gs = GeneSet(
            my_gene_set.id + str(i+1),
            my_gene_set.name,
            my_gene_set.genes
        )
        gene_sets.append(gs)
    return gene_sets


@pytest.fixture
def my_gene_set_coll(my_gene_sets):
    db = GeneSetCollection(my_gene_sets)
    return db


def test_init(my_gene_set_coll, my_gene_sets):
    assert isinstance(my_gene_set_coll, GeneSetCollection)
    assert my_gene_set_coll.gene_sets == my_gene_sets
    assert isinstance(repr(my_gene_set_coll), str)
    assert isinstance(str(my_gene_set_coll), str)
    assert isinstance(text(my_gene_set_coll), text)
    assert isinstance(my_gene_set_coll.hash, text)
    assert my_gene_set_coll.n == len(my_gene_sets)


# @pytest.mark.xfail(raises=ValueError, strict=True)  # not the right way
def test_id_clash(my_gene_set):
    """Make sure that GeneSetCollection complains about gene set ID clashes."""
    copy_ = deepcopy(my_gene_set)
    with pytest.raises(ValueError):
        db = GeneSetCollection([my_gene_set, copy_])


def test_lookup(my_gene_set_coll, my_gene_sets):
    for i, gs in enumerate(my_gene_sets):
        assert my_gene_set_coll[i] == gs
        assert my_gene_set_coll[gs.id] == gs
        assert my_gene_set_coll.get_by_index(i) == gs
        assert my_gene_set_coll.get_by_id(gs.id) == gs
        assert my_gene_set_coll.index(gs.id) == i


def test_tsv(my_gene_set_coll, tmpdir):
    tmp_file = str(tmpdir.join('gene_set.tsv'))
    my_gene_set_coll.write_tsv(tmp_file)
    other = GeneSetCollection.read_tsv(tmp_file)
    assert my_gene_set_coll == other