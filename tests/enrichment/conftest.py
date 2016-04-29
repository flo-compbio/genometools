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
from builtins import str as text

from string import ascii_lowercase

import pytest
import numpy as np
from xlmhg import xlmhg_test

from genometools.basic import GeneSet, GeneSetDB
from genometools.expression import ExpGenome
from genometools.enrichment import GSEAnalysis, GSEResult


@pytest.fixture
def my_genome():
    genes = [text(c) for c in ascii_lowercase]
    genome = ExpGenome.from_gene_names(genes)
    return genome


@pytest.fixture
def my_v():
    v = np.uint8([1, 0, 1, 1, 0, 1] + [0] * 12 + [1, 0])
    return v


@pytest.fixture
def my_ranked_genes(my_genome, my_v):
    return my_genome.genes[:my_v.size]


@pytest.fixture
def my_gene_set(my_ranked_genes, my_v):
    """Select the genes corresponding to the 1's in ``my_v``."""
    genes = [my_ranked_genes[i] for i in np.nonzero(my_v)[0]]
    gene_set = GeneSet('TestID', 'TestName', genes)
    return gene_set


@pytest.fixture
def my_uninteresting_gene_set(my_ranked_genes):
    """Select the last five ranked genes"""
    genes = my_ranked_genes[-5:]
    gene_set = GeneSet('BoringID', 'boring gene set', genes)
    return gene_set


@pytest.fixture
def my_gene_set_db(my_gene_set, my_uninteresting_gene_set):
    db = GeneSetDB([my_gene_set, my_uninteresting_gene_set])
    return db


@pytest.fixture
def my_result(my_v, my_ranked_genes, my_gene_set):
    N = my_v.size
    X = 1
    L = N
    sel = np.nonzero(my_v)[0]
    sel_genes = [my_ranked_genes[i] for i in sel]
    assert set(sel_genes) == my_gene_set.genes
    stat, n_star, pval = xlmhg_test(my_v, X, L)
    result = GSEResult(stat, n_star, pval, N, X, L, my_gene_set,
                       sel, sel_genes)
    return result

    # @pytest.fixture
# def my_result(my_v, my_gene_set):
#    N = my_v.size
#    stat, n_star, pval = xlmhg_test(my_v, X=1, L=N)
#    GSEResult(stat, n_star, pval, N, X, L, my_gene_set)