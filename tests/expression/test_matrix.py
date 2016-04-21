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
from builtins import *

import hashlib

import pytest
import numpy as np

from genometools.expression import ExpMatrix
from genometools.expression import ExpProfile


@pytest.fixture
def E_test():
    genes = ['a', 'b', 'c', 'd']
    samples = ['s1', 's2', 's3']
    X = np.arange(12, dtype=np.float64).reshape(4, 3)
    E = ExpMatrix(genes=genes, samples=samples, X=X)
    return E


def test_slice(E_test):
    e = E_test.iloc[:, 0]
    assert isinstance(e, ExpProfile)


def test_write_read(tmpdir, E_test):
    path = tmpdir.join('test.txt')
    E_test.write_tsv(str(path))
    data = open(str(path), mode='rb').read()
    h = hashlib.md5(data).hexdigest()
    assert h == 'd34bf3d376eb613e4fea894f7c9d601f'
    E2 = ExpMatrix.read_tsv(str(path))
    assert isinstance(E2, ExpMatrix)
    assert E_test.index.equals(E2.index)
    assert E_test.columns.equals(E2.columns)
    assert E_test.equals(E2)
