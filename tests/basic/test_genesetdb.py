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

"""Tests for the `GeneSetDB` class."""

# currently missing: test of the `GeneSetDB.read_msigdb_xml` class method

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)

from copy import deepcopy

import pytest
from genometools.basic import GeneSet, GeneSetDB


@pytest.fixture
def my_gene_sets(my_gene_set):
    gene_sets = []
    for i in range(3):
        gs = deepcopy(my_gene_set)
        gs.id += str(i+1)
        gene_sets.append(gs)
    return gene_sets


@pytest.fixture
def my_gene_set_db(my_gene_sets):
    db = GeneSetDB(my_gene_sets)
    return db


def test_init(my_gene_set_db, my_gene_sets):
    assert isinstance(my_gene_set_db, GeneSetDB)
    assert my_gene_set_db.gene_sets == my_gene_sets
    assert isinstance(repr(my_gene_set_db), str)
    assert isinstance(str(my_gene_set_db), str)
    assert isinstance(my_gene_set_db.hash, str)
    assert my_gene_set_db.n == len(my_gene_sets)


# @pytest.mark.xfail(raises=ValueError, strict=True)  # not the right way
def test_id_clash(my_gene_set):
    """Make sure that GeneSetDB complains about gene set ID clashes."""
    copy_ = deepcopy(my_gene_set)
    with pytest.raises(ValueError):
        db = GeneSetDB([my_gene_set, copy_])


def test_lookup(my_gene_set_db, my_gene_sets):
    for i, gs in enumerate(my_gene_sets):
        assert my_gene_set_db[i] == gs
        assert my_gene_set_db[gs.id] == gs
        assert my_gene_set_db.get_by_index(i) == gs
        assert my_gene_set_db.get_by_id(gs.id) == gs
        assert my_gene_set_db.index(gs.id) == i


def test_tsv(my_gene_set_db, tmpdir):
    tmp_file = str(tmpdir.join('gene_set.tsv'))
    my_gene_set_db.write_tsv(tmp_file)
    other = GeneSetDB.read_tsv(tmp_file)
    assert my_gene_set_db == other