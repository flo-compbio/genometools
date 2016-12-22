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

"""Utility functions for expression visualizations."""

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import *

import numpy as np


def read_colorscale(cmap_file):
    """Return a colorscale in the format expected by plotly.

    Parameters
    ----------
    cmap_file : str
        Path of a plain-text file containing the colorscale. 

    Returns
    -------
    list
        The colorscale.
        
    Notes
    -----
    A plotly colorscale is a list where each item is a pair
    (i.e., a list with two elements) consisting of
    a decimal number x between 0 and 1 and a corresponding "rgb(r,g,b)" string,
    where r, g, and b are integers between 0 and 255.

    The `cmap_file` is a tab-separated text file containing four columns
    (x,r,g,b), so that each row corresponds to an entry in the list
    described above.
    """
    assert isinstance(cmap_file, str)

    cm = np.loadtxt(cmap_file, delimiter='\t', dtype=np.float64)
    # x = cm[:, 0]
    rgb = np.int64(cm[:, 1:])  # normalize to 0-1?
    n = cm.shape[0]
    colorscale = []
    for i in range(n):
        colorscale.append(
            [i / float(n - 1),
             'rgb(%d, %d, %d)' % (rgb[i, 0], rgb[i, 1], rgb[i, 2])]
        )
    return colorscale

