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

"""Tests for the `GeneSet` class."""

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import str as text


import copy

import pytest

from genometools.ontology import GOTerm, GeneOntology


@pytest.fixture
def my_term():
    term = GOTerm('GO:0000000', 'regulation of test process',
                  'biological_process', 'This is a test GO term.')
    return term


@pytest.fixture
def my_other_term():
    term = GOTerm('GO:0000001', 'positive regulation of test process',
                  'biological_process', 'This is another test GO term.',
                  is_a=['GO:0000000'])
    return term


def test_basic(my_term, my_other_term):

    ontology = GeneOntology([my_term, my_other_term])

    assert isinstance(ontology, GeneOntology)
    assert isinstance(repr(ontology), str)
    assert isinstance(str(ontology), str)
    assert isinstance(text(ontology), text)
    assert isinstance(ontology.hash, text)

    # test access methods
    assert len(ontology) == 2
    assert my_term.id in ontology
    assert ontology[my_term.id] == my_term
    del ontology[my_term.id]
    assert my_term.id not in ontology
    ontology[my_term.id] = my_term
    assert my_term.id in ontology

    # test additional access methods
    assert ontology.get_term_by_id(my_term.id) == my_term
    assert ontology.get_term_by_acc(my_term.acc) == my_term

    # test comparisons
    other = copy.deepcopy(ontology)
    assert other == ontology
    del other[my_other_term.id]
    assert other != ontology

    # test iteration
    assert set(list(iter(ontology))) == set([my_term, my_other_term])

def test_real(my_gene_ontology_file):
    ontology = GeneOntology.read_obo(my_gene_ontology_file)
    assert isinstance(ontology, GeneOntology)