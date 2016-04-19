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

"""Methods for filtering expression data."""

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import *

import logging

import numpy as np

from . import ExpMatrix

logger = logging.getLogger(__name__)

def filter_variance(E, top):
    assert isinstance(E, ExpMatrix)
    assert isinstance(top, int)

    if top >= E.p:
        logger.warning('Variance filter with `top` parameter that is >= '
                       'the number of genes!')
        top = E.p

    a = np.argsort(np.var(E.X, axis=1))
    a = a[::-1]
    
    sel = np.zeros(E.p, dtype=np.bool_)
    sel[a[:top]] = True
    sel = np.nonzero(sel)[0]

    E = E.iloc[sel]
    return E
