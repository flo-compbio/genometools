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

import pytest
import numpy as np

from genometools.expression import ExpGene, ExpGenome, ExpMatrix


@pytest.fixture
def my_gene_names():
    return ['d', 'c', 'b', 'a']


@pytest.fixture
def my_unknown_gene_name():
    return 'X'


@pytest.fixture
def my_genes(my_gene_names):
    chromosomes = ['1', 'X', None, None]
    ensembl_ids = [None, 'ENSG0000001', 'ENSG0000002', None]
    return [ExpGene(g, chromosome=c, ensembl_id=e)
            for g, c, e in zip(my_gene_names, chromosomes, ensembl_ids)]


@pytest.fixture
def my_genome(my_genes):
    return ExpGenome(my_genes)


@pytest.fixture
def my_samples():
    return ['s1', 's2', 's3']


@pytest.fixture
def my_x():
    a = np.arange(4, dtype=np.float64)
    return a


@pytest.fixture
def my_X(my_x):
    X = []
    for i in range(0, -3, -1):
        X.append(np.roll(my_x,i))
    X = np.float64(X).T
    return X


@pytest.fixture
def my_matrix(my_gene_names, my_samples, my_X):
    #genes = ['a', 'b', 'c', 'd']
    #samples = ['s1', 's2', 's3']
    # X = np.arange(12, dtype=np.float64).reshape(4, 3)
    matrix = ExpMatrix(genes=my_gene_names, samples=my_samples, X=my_X)
    return matrix