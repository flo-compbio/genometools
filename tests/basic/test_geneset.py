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
from __future__ import (absolute_import, division,
                        print_function, unicode_literals)

import pytest

from genometools.basic import GeneSet


@pytest.fixture
def my_string():
    return 'Hello World'


def test_init(my_gene_set, my_gene_set2, my_genes):
    print('Enter')
    for gs in [my_gene_set, my_gene_set2]:
        assert isinstance(gs, GeneSet)
        assert isinstance(repr(gs), str)
        assert isinstance(str(gs), str)
        assert gs.size == len(my_genes)
    assert my_gene_set != my_gene_set2


def test_list(my_gene_set):
    l = my_gene_set.to_list()
    assert isinstance(l, list)
    assert len(l) == 6
    other = GeneSet.from_list(l)
    assert other == my_gene_set