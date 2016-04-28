from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import *

import pytest
import numpy as np

from genometools.expression import ExpGene, ExpGenome


@pytest.fixture(scope='session')
def my_genes():
    return ['d', 'c', 'b', 'a']


@pytest.fixture(scope='session')
def not_my_gene():
    return 'X'


@pytest.fixture(scope='session')
def my_exp_genes(my_genes):
    chromosomes = [['1'], ['X', 'Y'], None, None]
    ensembl_ids = [None, ['ENSG0000001'], ['ENSG0000002'], None]
    return [ExpGene(g, chromosomes=c, ensembl_ids=e)
            for g, c, e in zip(my_genes, chromosomes, ensembl_ids)]


@pytest.fixture(scope='session')
def my_genome(my_exp_genes):
    return ExpGenome(my_exp_genes)


@pytest.fixture(scope='session')
def my_x():
    a = np.arange(4, dtype=np.float64)
    return a


@pytest.fixture(scope='session')
def my_samples():
    return ['s1', 's2', 's3']
