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

"""Tests for the `GSEResult` clsas."""

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import str as text

from copy import deepcopy

# import pytest
import numpy as np

from genometools.enrichment import StaticGSEResult, RankBasedGSEResult


def test_ranked_basic(my_rank_based_result, my_v):
    assert isinstance(my_rank_based_result, RankBasedGSEResult)
    assert isinstance(repr(my_rank_based_result), str)
    assert isinstance(str(my_rank_based_result), str)
    assert isinstance(text(my_rank_based_result), text)
    assert isinstance(my_rank_based_result.hash, text)

    other = deepcopy(my_rank_based_result)
    assert other is not my_rank_based_result
    assert other == my_rank_based_result
    other.ind_genes[0] += 'Hello'
    assert other != my_rank_based_result

    assert my_rank_based_result.K == np.nonzero(my_v)[0].size
    assert my_rank_based_result.k == \
            int(np.sum(np.nonzero(my_v)[0] < \
            my_rank_based_result.cutoff))


def test_ranked_format(my_rank_based_result):
    assert isinstance(my_rank_based_result.get_pretty_format(), text)


def test_static_basic(my_static_result, my_v, my_static_genes):
    assert isinstance(my_static_result, StaticGSEResult)
    assert isinstance(repr(my_static_result), str)
    assert isinstance(str(my_static_result), str)
    assert isinstance(text(my_static_result), text)
    assert isinstance(my_static_result.hash, text)

    other = deepcopy(my_static_result)
    assert other is not my_static_result
    assert other == my_static_result
    other.N += 1
    assert other != my_static_result

    assert my_static_result.K == np.nonzero(my_v)[0].size
    assert my_static_result.n == len(my_static_genes)
    assert my_static_result.k == len(my_static_result.selected_genes)



def test_static_format(my_static_result):
    assert isinstance(my_static_result.get_pretty_format(), text)
