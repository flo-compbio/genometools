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
def e():
    genes = ['a','b','c','d']
    label = 's1'
    x = np.arange(4, dtype=np.float64)
    e = ExpProfile(genes=genes, label=label, x=x)
    return e

def test_expanddim(e):
    E = e.to_frame()
    assert isinstance(E, ExpMatrix)

def test_write_read(tmpdir, e):
    path = tmpdir.join('test.txt')
    e.write_tsv(str(path))
    data = open(str(path), mode='rb').read()
    h = hashlib.md5(data).hexdigest()
    assert h == '0edbe33c2c35354019a71cd85f11137d'
    e2 = ExpProfile.read_tsv(str(path))
    assert isinstance(e2, ExpProfile)
    assert e.name == e2.name
    assert e.index.equals(e2.index)
    assert e.equals(e2)
