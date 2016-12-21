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

"""Module containing the `SampleCorrelationHeatmap` class.

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
from . import read_colorscale

# from . import HeatmapSampleAnnotation, HeatmapBlockAnnotation

from collections import Iterable

logger = logging.getLogger(__name__)


class SampleCorrelationHeatmap(object):
    """A sample correlation heatmap."""

    _default_cmap_file = genometools._root.rstrip(os.sep) + os.sep + \
                         os.sep.join(['data', 'RdBu_r_colormap.tsv'])

    def __init__(self, corr_matrix,
                 sample_annotations=None,
                 block_annotations=None,
                 colorscale=None, colorbar_label=None):

        if sample_annotations is None:
            sample_annotations = []

        if block_annotations is None:
            block_annotations = []

        if colorscale is None:
            # use default colorscale
            colorscale = read_colorscale(self._default_cmap_file)

        assert isinstance(corr_matrix, ExpMatrix)
        assert isinstance(sample_annotations, Iterable)
        assert isinstance(block_annotations, Iterable)
        assert isinstance(colorscale, Iterable)
        if colorbar_label is not None:
            assert isinstance(colorbar_label, str)

        # make sure correlation matrix is square
        assert corr_matrix.shape[0] == corr_matrix.shape[1]

        self.corr_matrix = corr_matrix
        self.sample_annotations = sample_annotations
        self.block_annotations = block_annotations
        self.colorscale = colorscale
        self.colorbar_label = colorbar_label

    def get_figure(self, **kwargs):
        """Get a plotly figure of the heatmap."""

        emin = kwargs.pop('emin', -1.0)
        emax = kwargs.pop('emax', 1.0)
        width = kwargs.pop('width', 800)
        height = kwargs.pop('height', 600)
        margin_left = kwargs.pop('margin_left', 100)
        margin_bottom = kwargs.pop('margin_bottom', 60)
        margin_top = kwargs.pop('margin_top', 30)
        margin_right = kwargs.pop('margin_right', 0)
        colorbar_size = kwargs.pop('colorbar_size', 0.4)

        xaxis_label = kwargs.pop('xaxis_label', None)
        yaxis_label = kwargs.pop('yaxis_label', None)
        xticks_angle = kwargs.pop('xaxis_angle', 30)
        font = kwargs.pop('font', '"Droid Serif", "Open Serif", serif')
        font_size = kwargs.pop('font_size', 12)
        title = kwargs.pop('title', None)
        title_font_size = kwargs.pop('title_font_size', None)
        annotation_font_size = kwargs.pop('annotation_font_size', None)
        show_sample_labels = kwargs.pop('show_sample_labels', 'x')

        if show_sample_labels not in ['none', 'x', 'y', 'both']:
            raise ValueError('"show_sample_labels" must be "none", "x", "y", '
                             'or "both".')

        padding_top = kwargs.pop('padding_top', 0.1)
        padding_right = kwargs.pop('padding_top', 0.1)

        xaxis_nticks = kwargs.pop('xaxis_nticks', None)
        #yaxis_nticks = kwargs.pop('yaxis_nticks', None)

        if title_font_size is None:
            title_font_size = font_size

        if annotation_font_size is None:
            annotation_font_size = font_size

        colorbar_label = self.colorbar_label or 'Pearson Correlation'

        ### set up heatmap

        colorbar = go.ColorBar(
            lenmode='fraction',
            len=colorbar_size,
            title=colorbar_label,
            titlefont = dict(
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
                except ValueError:
                    fixed_labels.append(str(l))
                else:
                    fixed_labels.append(str(l) + '_')
            return fixed_labels

        x = fix_plotly_label_bug(self.corr_matrix.samples)
        y = x

        data = [
            go.Heatmap(
                z=self.corr_matrix.X,
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

        xshowticklabels = False
        yshowticklabels = False

        ### set up layout
        if show_sample_labels == 'x':
            xshowticklabels = True
        elif show_sample_labels == 'y':
            yshowticklabels = True
        elif show_sample_labels == 'both':
            xshowticklabels = True
            yshowticklabels = True

        xticks = 'outside'
        yticks = 'outside'

        if xaxis_label is None:
            if self.corr_matrix.samples.name is not None:
                xaxis_label = self.corr_matrix.samples.name
            else:
                xaxis_label = 'Samples'
            xaxis_label = xaxis_label + ' (n = %d)' % self.corr_matrix.n

        if yaxis_label is None:
            yaxis_label = xaxis_label

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
                showticklabels=xshowticklabels,
                ticks=xticks,
                nticks=xaxis_nticks,
                tickangle=xticks_angle,
                #range=[-0.5, self.corr_matrix.n-0.5],
                showline=True,
                zeroline=False,
                showgrid=False,
            ),
            yaxis=go.YAxis(
                title=yaxis_label,
                titlefont=dict(size=title_font_size),
                showticklabels=yshowticklabels,
                ticks=xticks,
                nticks=xaxis_nticks,
                autorange='reversed',
                showline=True,
                zeroline=False,
                showgrid=False,
            ),

            margin=go.Margin(
                l=margin_left,
                t=margin_top,
                b=margin_bottom,
                r=margin_right,
                pad=0
            ),
        )

        ### add annotations

        # we need separate, but overlaying, axes to place the annotations
        layout['xaxis2'] = go.XAxis(
            overlaying='x',
            showline=False,
            tickfont=dict(size=0),
            autorange=False,
            #range=[-0.5, self.corr_matrix.n-0.5],
            range=[-0.5, self.corr_matrix.n-0.5],
            ticks='',
            showticklabels=False,
            zeroline=False,
            showgrid=False,
        )

        layout['yaxis2'] = go.YAxis(
            overlaying='y',
            showline=False,
            tickfont=dict(size=0),
            autorange=False,
            range=[self.corr_matrix.n-0.5, -0.5],
            ticks='',
            showticklabels=False,
            zeroline=False,
            showgrid=False,
        )

        # generate coordinates and labels for the block annotations
        k = len(self.block_annotations)
        block_coords = np.zeros((k, 2), dtype=np.float64)
        block_labels = []
        for i, ann in enumerate(self.block_annotations):
            block_coords[i, :] = [ann.start_index-0.5, ann.end_index+0.5]
            block_labels.append(ann.label)

        # this produces the squares for the block annotations
        for i in range(k):
            mn = block_coords[i, 0]
            mx = block_coords[i, 1]
            data.append(
                go.Scatter(
                    x=[mn, mx, mx, mn, mn],
                    y=[mn, mn, mx, mx, mn],
                    mode='lines',
                    hoverinfo='none',
                    showlegend=False,
                    line=dict(color='black'),
                    xaxis='x2',
                    yaxis='y2',
                )
            )

        # - this produces the square labels for the block annotations
        # - we use plotly annotations, so that the labels are not limited
        #   by the plotting area
        for i in range(k):
            mn = block_coords[i, 0]
            mx = block_coords[i, 1]
            layout.annotations.append(
                dict(
                    x=mx,
                    y=(mn+mx)/2.0,
                    text=block_labels[i],
                    xref='x2',
                    yref='y2',
                    showarrow=False,
                    #ax=20,
                    #ay=-20,
                    xanchor='left',
                    yanchor='middle',
                    font=dict(
                        color='black',
                        size=annotation_font_size,
                        family='serif bold',
                    ),
                    bgcolor='white',
                    opacity=0.7,
                    #textanchor='top right'
                )
            )

        fig = go.Figure(
            data=data,
            layout=layout
        )

        return fig
