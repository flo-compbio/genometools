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

import pytest
# import numpy as np
from plotly import graph_objs as go

from genometools.expression import cluster, ExpMatrix
from genometools.expression.visualize import *


@pytest.fixture
def my_gene_annotation():
    return HeatMapGeneAnnotation('b', 'blue', 'Important gene')


@pytest.fixture
def my_sample_annotation():
    return HeatMapSampleAnnotation('s2', 'green', 'Important sample')


@pytest.fixture
def my_heatmap(my_matrix, my_gene_annotation, my_sample_annotation):
    assert isinstance(my_matrix, ExpMatrix)
    assert isinstance(my_gene_annotation, HeatMapGeneAnnotation)
    assert isinstance(my_sample_annotation, HeatMapSampleAnnotation)
    heatmap = ExpHeatMap(
        my_matrix,
        gene_annotations=[my_gene_annotation],
        sample_annotations=[my_sample_annotation]
    )
    return heatmap


def test_basic(my_heatmap):
    assert isinstance(my_heatmap, ExpHeatMap)
    # assert isinstance(repr(my_heatmap), str)
    # assert isinstance(str(my_heatmap), str)
    # assert isinstance(text(my_heatmap), text)
    # assert isinstance(my_heatmap.hash, text)

def test_figure(my_heatmap):
    figure = my_heatmap.get_figure()
    assert isinstance(figure, go.Figure)