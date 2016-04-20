# Copyright (c) 2016 Florian Wagner

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

"""Methods for plotting expression data."""

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import *

import os
import logging

import pandas as pd
import numpy as np

import plotly.graph_objs as go

from .. import _root
from .. import misc
from . import ExpMatrix

logger = logging.getLogger(__name__)

default_cmap_file = _root.rstrip(os.sep) + os.sep + \
                    os.sep.join(['data', 'RdBu_r_colormap.tsv'])

def _read_colorscale(cmap_file):
    """Return a colorscale in the format that plotly expects it.

    Specifically, the scale should be a list containing pairs consisting of
    a normalized value x (between 0 and 1) and a corresponding "rgb(r,g,b)"
    string, where r,g,b are integers from 0 to 255.

    The ``cmap_file`` is a tab-separated text file containing four columns
    (x,r,g,b), so that each row corresponds to an entry in the list described
    above.
    """
    assert isinstance(cmap_file, str)

    cm = np.loadtxt(cmap_file, delimiter = '\t', dtype = np.float64)
    x = cm[:,0]
    rgb = np.int64(cm[:,1:]) # normalize to 0-1?
    n = cm.shape[0]
    colorscale = []
    for i in range(n):
        colorscale.append(
            [i / float(n - 1),
             'rgb(%d, %d, %d)' %(rgb[i,0], rgb[i,1], rgb[i,2])]
        )
    return colorscale

def get_heatmap(
        E, colorscale=None, title=None, emin=None, emax=None,
        width=800, height=400,
        margin_left=100, margin_bottom=60, margin_top=30,
        colorbar_label='Expression', colorbar_size=0.4,
        yaxis_label = 'Genes', xaxis_label = 'Samples',
        yaxis_nticks = None, xaxis_nticks = None,
        xtick_angle = 30,
        font = '"Droid Serif", "Open Serif", serif',
        font_size = 12, title_font_size = None,
        show_sample_labels = True
    ):
    
    assert isinstance(E, ExpMatrix)

    if colorscale is None:
        # load default colorscale
        colorscale = _read_colorscale(default_cmap_file)

    # emin and/or emax are unspecified, set to data min/max values
    if emax is None:
        emax = E.X.max()
    if emin is None:
        emin = E.X.min()

    if title_font_size is None:
        title_font_size = font_size

    colorbar = dict(
        lenmode = 'fraction',
        len = colorbar_size,
        title = colorbar_label,
        titleside = 'right',
        xpad = 0,
        ypad = 0,
        outlinewidth = 0, # no border
        thickness = 20, # in pixels
        #outlinecolor = '#000000',
    )

    data = [
        go.Heatmap(
            z = E.X,
            x = E.samples,
            y = E.genes,
            zmin = emin,
            zmax = emax,
            colorscale = colorscale,
            colorbar = colorbar,
            hoverinfo = 'x+y+z',
        ),
    ]

    xticks = 'outside'
    if not show_sample_labels:
        xticks = ''

    if xaxis_label is not None:
        xaxis_label = xaxis_label + ' (n = %d)' %(E.n)

    layout = go.Layout(
        width = width,
        height = height,
        title = title,
        titlefont = dict(size = title_font_size),

        font = dict(size = font_size, family = font),

        xaxis = dict(title = xaxis_label,
                     titlefont = dict(size = title_font_size),
                     showticklabels = show_sample_labels, ticks = xticks,
                     nticks = xaxis_nticks, tickangle = xtick_angle,
                     showline = True),

        yaxis = dict(title = yaxis_label,
                     titlefont = dict(size = title_font_size),
                     nticks = yaxis_nticks, autorange = 'reversed',
                     showline = True),

        margin = dict(l = margin_left, t = margin_top, b = margin_bottom, r = 0,
                      pad = 0),
    )

    fig = go.Figure(data = data, layout = layout)

    return fig
