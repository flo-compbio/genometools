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

from genometools.enrichment import GSEResult


def test_init(my_result, my_v):
    assert isinstance(my_result, GSEResult)
    assert isinstance(repr(my_result), str)
    assert isinstance(str(my_result), str)
    assert isinstance(text(my_result), text)
    assert isinstance(my_result.hash, text)

    other = deepcopy(my_result)
    assert other is not my_result
    assert other == my_result
    other.ind_genes[0] += 'Hello'
    assert other != my_result

    assert my_result.K == np.nonzero(my_v)[0].size
    assert my_result.k == int(np.sum(np.nonzero(my_v)[0] < my_result.cutoff))


def test_format(my_result):
    assert isinstance(my_result.get_pretty_format(), text)