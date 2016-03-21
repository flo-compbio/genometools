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

import logging

import pandas as pd
import numpy as np

from bokeh.models import HoverTool

from .. import misc
from . import ExpMatrix

logger = logging.getLogger(__name__)

def plot_heatmap(E, cmap, title = None, vmin = None, vmax = None, width = 800, height = 400, yaxis_label = 'Genes', xaxis_label = 'Samples',
                 font = 'Computer Modern Roman', show_sample_labels = False, font_size = '10pt'):
    from bokeh.plotting import ColumnDataSource, figure, show
    
    # create color mapping
    if vmax is None:
        vmax = E.X.max()
    if vmin is None:
        vmin = E.X.min()
    C = np.int64(np.round(255 * (E.X - vmin) / (vmax - vmin)))
    C[C<0] = 0
    C[C>255] = 255
    
    genes = [g.replace(':','.') for g in tuple(E.genes)]
    samples = [s.replace(':','.') for s in tuple(E.samples)]
    
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
    hover = p.select(dict(type=HoverTool))
    hover.tooltips = [
        (yaxis_label[:-1] + ' (y)', '$y'),
        (xaxis_label[:-1] + ' (x)', '$x'),
        ('Value', '@data_expr'),
    ]
    
    hm = show(p)
    return hm
