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

"""Tests for functions in the `filter` module."""

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import str as text

# import numpy as np

from genometools.expression import filter


def test_init(my_matrix):
    #print(np.mean(my_matrix.X, axis=1))
    #print(np.var(my_matrix.X, ddof=1, axis=1))
    #print(np.percentile(my_matrix.X, 25, axis=1))
    top = 2
    percentile = 25
    other = filter.filter_variance(my_matrix, top)
    assert other.shape[0] == top
    other = filter.filter_mean(my_matrix, top)
    assert other.shape[0] == top
    other = filter.filter_percentile(my_matrix, top, percentile=percentile)
    assert other.shape[0] == top