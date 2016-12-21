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

"""Module containing the `ExpHeatmap` class.

"""

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
_oldstr = str
from builtins import *

import os
import logging

from plotly import graph_objs as go
import numpy as np

import genometools
from .. import ExpMatrix
from . import read_colorscale

from collections import Iterable

logger = logging.getLogger(__name__)


class ExpHeatmap(object):
    """An expression heatmap."""

    # TODO: docstrings, __str__, __repr__, hash

    _default_cmap_file = genometools._root.rstrip(os.sep) + os.sep + \
                         os.sep.join(['data', 'RdBu_r_colormap.tsv'])

    def __init__(self, matrix,
                 gene_annotations=None, sample_annotations=None,
                 colorscale=None, colorbar_label=None, title=None):

        if gene_annotations is None:
            gene_annotations = []

        if sample_annotations is None:
            sample_annotations = []

        if colorscale is None:
            # use default colorscale
            colorscale = read_colorscale(self._default_cmap_file)

        assert isinstance(matrix, ExpMatrix)
        assert isinstance(gene_annotations, Iterable)
        assert isinstance(sample_annotations, Iterable)
        assert isinstance(colorscale, Iterable)
        if colorbar_label is not None:
            assert isinstance(colorbar_label, (str, _oldstr))

        self.matrix = matrix
        self.gene_annotations = gene_annotations
        self.sample_annotations = sample_annotations
        self.colorscale = colorscale
        self.colorbar_label = colorbar_label
        self.title = title

    def get_figure(
            self, title=None, emin=None, emax=None,
            width=800, height=400,
            margin_left=100, margin_bottom=60, margin_top=30, margin_right=0,
            colorbar_size=0.4,
            xaxis_label=None, yaxis_label=None,
            xaxis_nticks=None, yaxis_nticks=None,
            xtick_angle=30,
            font='"Droid Serif", "Open Serif", serif',
            font_size=12, title_font_size=None,
            show_sample_labels=True, **kwargs):
        """Generate a plotly figure of the heatmap."""

        # emin and/or emax are unspecified, set to data min/max values
        if emax is None:
            emax = self.matrix.X.max()
        if emin is None:
            emin = self.matrix.X.min()

        if title is None:
            title = self.title

        if title_font_size is None:
            title_font_size = font_size

        colorbar_label = self.colorbar_label or 'Expression'

        colorbar = go.ColorBar(
            lenmode='fraction',
            len=colorbar_size,
            title=colorbar_label,
            titlefont=dict(
                size=title_font_size,
            ),
            titleside='right',
            xpad=0,
            ypad=0,
            outlinewidth=0,  # no border
            thickness=20,  # in pixels
            # outlinecolor = '#000000',
        )

        def fix_plotly_label_bug(labels):
            """
            This fixes a bug whereby plotly treats labels that look
            like numbers (integers or floats) as numeric instead of
            categorical, even when they are passed as strings. The fix consists
            of appending an underscore to any label that looks like a number.
            """
            assert isinstance(labels, Iterable)
            fixed_labels = []
            for l in labels:
                try:
                    float(l)
                except (ValueError, TypeError):
                    fixed_labels.append(str(l))
                else:
                    fixed_labels.append(str(l) + '_')
            return fixed_labels

        x = fix_plotly_label_bug(self.matrix.samples)
        y = fix_plotly_label_bug(self.matrix.genes)

        data = [
            go.Heatmap(
                z=self.matrix.X,
                x=x,
                y=y,
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

        if xaxis_label is None:
            if self.matrix.samples.name is not None:
                xaxis_label = self.matrix.samples.name
            else:
                xaxis_label = 'Samples'
            xaxis_label = xaxis_label + ' (n = %d)' % self.matrix.n

        if yaxis_label is None:
            if self.matrix.genes.name is not None:
                yaxis_label = self.matrix.genes.name
            else:
                yaxis_label = 'Genes'
            yaxis_label = yaxis_label + ' (p = %d)' % self.matrix.p

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
                r=margin_right,
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
                    x=[xmn, xmx, xmx, xmn, xmn],
                    y=[ymn, ymn, ymx, ymx, ymn],
                    mode='lines',
                    hoverinfo='none',
                    showlegend=False,
                    line=dict(color=ann.color),
                    xaxis='x2',
                    yaxis='y2',
                    opacity=0.5,
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
                        opacity=1-ann.transparency,
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
                    x=[xmn, xmx, xmx, xmn, xmn],
                    y=[ymn, ymn, ymx, ymx, ymn],
                    mode='lines',
                    hoverinfo='none',
                    showlegend=False,
                    line=dict(color=ann.color),
                    xaxis='x2',
                    yaxis='y2',
                    opacity=1.0)
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
                        opacity=1-ann.transparency,
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
