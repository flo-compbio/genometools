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

"""Tests for the `GSEAnalysis` class."""

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import str as text

from string import ascii_lowercase

import pytest

# from genometools.expression import ExpGenome
from genometools import misc
from genometools.enrichment import GeneSetEnrichmentAnalysis

logger = misc.get_logger('genometools', verbose=True)


@pytest.fixture
def my_analysis(my_genome, my_gene_set_coll):
    analysis = GeneSetEnrichmentAnalysis(my_genome, my_gene_set_coll)
    return analysis


def test_basic(my_analysis, my_genome):
    assert isinstance(my_analysis, GeneSetEnrichmentAnalysis)
    assert isinstance(repr(my_analysis), str)
    assert isinstance(str(my_analysis), str)
    assert isinstance(text(my_analysis), text)

    assert my_analysis.genome is not my_genome
    assert len(my_analysis.genome) == len(my_genome)


def test_rank_based_analysis(my_analysis, my_ranked_genes,
                             my_uninteresting_gene_set):
    """Tests rank-based gene set enrichment analysis."""
    pval_thresh = 0.025
    X_frac = 0
    X_min = 1
    L = len(my_ranked_genes)

    # test if rank-based enrichment works
    enriched = my_analysis.get_rank_based_enrichment(
        my_ranked_genes, pval_thresh, X_frac, X_min, L,
        adjust_pval_thresh=False)
    assert isinstance(enriched, list)
    assert len(enriched) == 1

    # test if selective testing of individual gene sets works
    enriched = my_analysis.get_rank_based_enrichment(
        my_ranked_genes, pval_thresh, X_frac, X_min, L,
        adjust_pval_thresh=False,
        gene_set_ids=[my_uninteresting_gene_set.id])
    assert isinstance(enriched, list)
    assert len(enriched) == 0


def test_static_analysis(my_analysis, my_static_genes):
    """Tests static gene set enrichment analysis."""
    pval_thresh = 0.05

    enriched = my_analysis.get_static_enrichment(
        my_static_genes, pval_thresh, adjust_pval_thresh=False)
    assert isinstance(enriched, list)
    assert len(enriched) == 1