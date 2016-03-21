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

"""Methods for plotting expression data."""

import os
import logging

import pandas as pd
import numpy as np

from bokeh.plotting import ColumnDataSource, figure, show
from bokeh.models import HoverTool

from .. import _root
from .. import misc
from . import ExpMatrix

logger = logging.getLogger(__name__)

default_cmap_file = _root.rstrip(os.sep) + os.sep + \
                    os.sep.join(['data', 'RdBu_r_colormap.tsv'])

def plot_heatmap(E, cmap = None, title = None, vmin = None, vmax = None,
                 width = 800, height = 400,
                 yaxis_label = 'Genes', xaxis_label = 'Samples',
                 font = None, font_size = None, title_font_size = None,
                 show_sample_labels = None):
    
    if cmap is None:
        # load default colormap
        cmap = np.loadtxt(default_cmap_file)

    # vmin and/or vmax are unspecified, set to data min/max values
    if vmax is None:
        vmax = E.X.max()
    if vmin is None:
        vmin = E.X.min()

    if font is None:
        font = 'Computer Modern Roman'

    if font_size is None:
        font_size = '10pt'

    if title_font_size is None:
        title_font_size = '14pt'

    if show_sample_labels is None:
        show_sample_labels = False
    
    # create color mapping
    C = np.int64(np.round(255 * (E.X - vmin) / (vmax - vmin)))
    C[C<0] = 0
    C[C>255] = 255

    # colons are not allowed in category names
    genes = [g.replace(':', '.') for g in E.genes]
    samples = [s.replace(':', '.') for s in E.samples]
    
    data_genes = []
    data_samples = []
    data_colors = []
    data_expr = []
    for i, g in enumerate(genes):
        for j, s in enumerate(samples):
            data_genes.append(g)
            data_samples.append(s)
            col = cmap[C[i,j],:]
            data_colors.append('#%02x%02x%02x' % (col[0], col[1], col[2]))
            data_expr.append(E.X[i,j])
    
    tools = 'reset,wheel_zoom,pan,save,hover,box_zoom,undo,redo'
    p=figure(
        title = title,
        x_range = samples,
        y_range = list(reversed(genes)), # reverse Y axis
        tools = tools
    )
    
    #p.add_tools(HoverTool())
    
    source = ColumnDataSource(
        data = dict(
            data_genes = data_genes,
            data_samples = data_samples,
            data_colors = data_colors,
            data_expr = data_expr
        )
    )
    
    p.rect('data_samples', 'data_genes', 1.0, 1.0, source = source, color = data_colors, line_color = None,
          dilate = True)
    
    p.plot_width = width
    p.plot_height = height
    
    p.toolbar_location = 'left'
    
    p.grid.grid_line_color = None
    p.axis.axis_line_color = None
    
    p.title_text_font = font
    p.title_text_font_size = title_font_size
    p.axis.axis_label_text_font = font
    p.axis.major_label_text_font = font
    
    p.axis.major_label_text_font_size = font_size
    if not show_sample_labels:
        p.xaxis.major_tick_line_color = None
        p.xaxis.major_label_text_font_size = '0pt'

    p.xaxis.axis_label = xaxis_label + ' (n = %d)' %(len(samples))
    p.yaxis.axis_label = yaxis_label
    #p.yaxis.major_label_standoff = 100
    p.xaxis.major_label_orientation = np.pi/3
    
    p.min_border_right = 10
    p.min_border_top = 10
    p.min_border_bottom = 10
    
    # set hover tooltips
    hover = p.select(dict(type = HoverTool))
    hover.tooltips = [
        (xaxis_label[:-1] + ' (x)', '$x'),
        (yaxis_label[:-1] + ' (y)', '$y'),
        ('Value', '@data_expr'),
    ]
    
    hm = show(p)
    return hm
