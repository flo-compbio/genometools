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

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)

from genometools.basic import GeneSet

import pytest


@pytest.fixture
def my_genes():
    return ['a', 'b', 'c']


@pytest.fixture
def my_gene_set(my_genes):
    gene_set = GeneSet('TestID', 'TestName', my_genes,
                       source='TestSource', collection='TestCollection',
                       description='Test GeneSet.')
    return gene_set


@pytest.fixture
def my_gene_set2(my_genes):
    # a gene set with all optional attributes set to None
    gene_set = GeneSet('TestID', 'TestName', my_genes)
    return gene_set