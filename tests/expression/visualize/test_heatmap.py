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

"""Tests for functions in `cluster` module."""

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import str as text, int as newint

import os

import pytest
# import numpy as np
from plotly import graph_objs as go

from genometools.expression import cluster, ExpMatrix
from genometools.expression.visualize import *


@pytest.fixture
def my_gene_annotation():
    return HeatmapGeneAnnotation('b', 'blue',
                                 label='Important gene',
                                 transparency=0.2)

@pytest.fixture
def my_sample_annotation():
    return HeatmapSampleAnnotation('s2', 'green')


@pytest.fixture
def my_block_annotation():
    return HeatmapBlockAnnotation(0, 1, label='My label')


@pytest.fixture
def my_heatmap_kw(my_gene_annotation, my_sample_annotation):
    assert isinstance(my_gene_annotation, HeatmapGeneAnnotation)
    assert isinstance(my_sample_annotation, HeatmapSampleAnnotation)
    heatmap_kw = dict(
        gene_annotations=[my_gene_annotation],
        sample_annotations=[my_sample_annotation],
        colorbar_label='My colorbar',
        title='Heatmap title',
    )
    return heatmap_kw


@pytest.fixture
def my_heatmap(my_matrix, my_heatmap_kw):
    assert isinstance(my_matrix, ExpMatrix)
    assert isinstance(my_heatmap_kw, dict)
    heatmap = ExpHeatmap(my_matrix, **my_heatmap_kw)
    return heatmap


@pytest.fixture
def my_sample_corr_heatmap(my_matrix, my_block_annotation):
    return None


def test_colorscale():
    assert os.path.isfile(ExpHeatmap._default_cmap_file)


def test_basic(my_heatmap):
    assert isinstance(my_heatmap, ExpHeatmap)
    # assert isinstance(repr(my_heatmap), str)
    # assert isinstance(str(my_heatmap), str)
    # assert isinstance(text(my_heatmap), text)
    # assert isinstance(my_heatmap.hash, text)

def test_figure(my_heatmap):
    figure = my_heatmap.get_figure()
    assert isinstance(figure, go.Figure)



def test_matrix_heatmap(my_matrix, my_gene_annotation, my_sample_annotation,
                        my_heatmap_kw):
    heatmap = my_matrix.get_heatmap(
        highlight_genes=[my_gene_annotation.gene],
        highlight_samples=[my_sample_annotation.sample],
        highlight_color=my_gene_annotation.color)
    assert isinstance(heatmap, ExpHeatmap)

    heatmap = my_matrix.get_heatmap(**my_heatmap_kw)
    assert isinstance(heatmap, ExpHeatmap)



def test_matrix_figure(my_matrix, my_heatmap_kw):
    fig = my_matrix.get_figure(heatmap_kw=my_heatmap_kw)
    assert isinstance(fig, go.graph_objs.Figure)

#def test_matrix_figure(my_matrix, my_heatmap_kw):
#    fig = my_matrix.get_figure()