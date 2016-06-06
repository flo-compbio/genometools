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

"""Module containing the `ExpHeatMap` class.

"""

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import *

import os
import logging

from plotly import graph_objs as go
import numpy as np

import genometools
from .. import ExpMatrix

from collections import Iterable

logger = logging.getLogger(__name__)


class ExpHeatMap(object):

    # TODO: docstrings, __str__, __repr__, hash

    _default_cmap_file = genometools._root.rstrip(os.sep) + os.sep + \
                         os.sep.join(['data', 'RdBu_r_colormap.tsv'])

    def __init__(self, matrix,
                 gene_annotations=None, sample_annotations=None,
                 colorscale=None):

        if gene_annotations is None:
            gene_annotations = []

        if sample_annotations is None:
            sample_annotations = []

        if colorscale is None:
            # use default colorscale
            colorscale = self._read_colorscale(self._default_cmap_file)

        assert isinstance(matrix, ExpMatrix)
        assert isinstance(gene_annotations, Iterable)
        assert isinstance(sample_annotations, Iterable)
        assert isinstance(colorscale, Iterable)

        self.matrix = matrix
        self.gene_annotations = gene_annotations
        self.sample_annotations = sample_annotations
        self.colorscale = colorscale

    @staticmethod
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

        cm = np.loadtxt(cmap_file, delimiter='\t', dtype=np.float64)
        # x = cm[:, 0]
        rgb = np.int64(cm[:, 1:])  # normalize to 0-1?
        n = cm.shape[0]
        colorscale = []
        for i in range(n):
            colorscale.append(
                [i / float(n-1),
                 'rgb(%d, %d, %d)' % (rgb[i, 0], rgb[i, 1], rgb[i, 2])]
            )
        return colorscale

    def get_figure(
            self, title=None, emin=None, emax=None,
            width=800, height=400,
            margin_left=100, margin_bottom=60, margin_top=30,
            colorbar_label='Express`ion', colorbar_size=0.4,
            xaxis_label='Samples', yaxis_label='Genes',
            xaxis_nticks=None, yaxis_nticks=None,
            xtick_angle=30,
            font='"Droid Serif", "Open Serif", serif',
            font_size=12, title_font_size=None,
            show_sample_labels=True, **kwargs):

        # emin and/or emax are unspecified, set to data min/max values
        if emax is None:
            emax = self.matrix.X.max()
        if emin is None:
            emin = self.matrix.X.min()

        if title_font_size is None:
            title_font_size = font_size

        colorbar = go.ColorBar(
            lenmode='fraction',
            len=colorbar_size,
            title=colorbar_label,
            titleside='right',
            xpad=0,
            ypad=0,
            outlinewidth=0,  # no border
            thickness=20,  # in pixels
            # outlinecolor = '#000000',
        )

        data = [
            go.Heatmap(
                z=self.matrix.X,
                x=self.matrix.samples.tolist(),
                y=self.matrix.genes.tolist(),
                zmin=emin,
                zmax=emax,
                colorscale=self.colorscale,
                colorbar=colorbar,
                hoverinfo='x+y+z',
                **kwargs
            ),
        ]

        xticks = 'outside'
        if not show_sample_labels:
            xticks = ''

        if xaxis_label is not None:
            xaxis_label = xaxis_label + ' (n = %d)' % self.matrix.n

        layout = go.Layout(
            width=width,
            height=height,
            title=title,
            titlefont=go.Font(
                size=title_font_size
            ),
            font=go.Font(
                size=font_size,
                family=font
            ),
            xaxis=go.XAxis(
                title=xaxis_label,
                titlefont=dict(size=title_font_size),
                showticklabels=show_sample_labels,
                ticks=xticks,
                nticks=xaxis_nticks,
                tickangle=xtick_angle,
                showline=True
            ),
            yaxis=go.YAxis(
                title=yaxis_label,
                titlefont=dict(size=title_font_size),
                nticks=yaxis_nticks,
                autorange='reversed',
                showline=True
            ),

            margin=go.Margin(
                l=margin_left,
                t=margin_top,
                b=margin_bottom,
                r=0,
                pad=0
            ),
        )

        # add annotations

        # we need separate, but overlaying, axes to place the annotations
        layout['xaxis2'] = go.XAxis(
            overlaying = 'x',
            showline = False,
            tickfont = dict(size=0),
            autorange=False,
            range=[-0.5, self.matrix.n-0.5],
            ticks='',
            showticklabels=False
        )

        layout['yaxis2'] = go.YAxis(
            overlaying='y',
            showline=False,
            tickfont=dict(size=0),
            autorange=False,
            range=[self.matrix.p-0.5, -0.5],
            ticks='',
            showticklabels=False
        )

        # gene (row) annotations
        for ann in self.gene_annotations:
            i = self.matrix.genes.get_loc(ann.gene)
            xmn = -0.5
            xmx = self.matrix.n-0.5
            ymn = i-0.5
            ymx = i+0.5
            data.append(
                go.Scatter(
                    x=[xmn,xmx,xmx,xmn,xmn],
                    y=[ymn,ymn,ymx,ymx,ymn],
                    mode='lines',
                    hoverinfo='none',
                    showlegend=False,
                    line=dict(color=ann.color),
                    xaxis='x2',
                    yaxis='y2',
                    opacity=0.5
                )
            )
            if ann.label is not None:
                layout.annotations.append(
                    go.Annotation(
                        text=ann.label,
                        x=0.01,
                        y=i-0.5,
                        #y=i+0.5,
                        xref='paper',
                        yref='y2',
                        xanchor='left',
                        yanchor='bottom',
                        showarrow=False,
                        bgcolor='white',
                        opacity=0.6,
                        borderpad=0,
                        #textangle=30,
                        font=dict(color=ann.color)
                    )
                )

        # sample (column) annotations
        for ann in self.sample_annotations:
            j = self.matrix.samples.get_loc(ann.sample)
            xmn = j-0.5
            xmx = j+0.5
            ymn = -0.5
            ymx = self.matrix.p-0.5
            data.append(
                go.Scatter(
                    x=[xmn,xmx,xmx,xmn,xmn],
                    y=[ymn,ymn,ymx,ymx,ymn],
                    mode='lines',
                    hoverinfo='none',
                    showlegend=False,
                    line=dict(color=ann.color),
                    xaxis='x2',
                    yaxis='y2',
                    opacity=0.5)
            )
            if ann.label is not None:
                layout.annotations.append(
                    go.Annotation(
                        text=ann.label,
                        y=0.99,
                        x=j+0.5,
                        #y=i+0.5,
                        xref='x2',
                        yref='paper',
                        xanchor='left',
                        yanchor='top',
                        showarrow=False,
                        bgcolor='white',
                        opacity=0.6,
                        borderpad=0,
                        textangle=90,
                        font=dict(color=ann.color)
                    )
                )

        fig = go.Figure(
            data=data,
            layout=layout
        )

        return fig