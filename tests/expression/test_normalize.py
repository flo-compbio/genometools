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

import numpy as np

from genometools.expression import ExpMatrix
from genometools.expression import normalize


def test_quantile_normalize(my_matrix):
    normal = normalize.quantile_normalize(my_matrix, inplace=False)
    assert isinstance(normal, ExpMatrix)
    assert normal.shape == my_matrix.shape
    assert np.all(normal.index == my_matrix.index)
    assert np.all(normal.columns == my_matrix.columns)

    X = np.sort(normal.X, axis=0)
    assert np.allclose(X, np.tile(np.mean(X, axis=1), (X.shape[1], 1)).T)