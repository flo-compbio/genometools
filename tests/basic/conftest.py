from __future__ import (absolute_import, division,
                        print_function, unicode_literals)

from genometools.basic import GeneSet

import pytest


@pytest.fixture(scope='session')
def my_genes():
    return ['a', 'b', 'c']


@pytest.fixture(scope='session')
def my_gene_set(my_genes):
    gene_set = GeneSet('TestID', 'TestName', my_genes, source='TestSource',
                 collection='TestCollection', description='Test GeneSet.')
    return gene_set


@pytest.fixture(scope='session')
def my_gene_set2(my_genes):
    # a gene set with all optional attributes set to None
    gene_set = GeneSet('TestID', 'TestName', my_genes)
    return gene_set