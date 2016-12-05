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
def my_other_term():
    go_term = GOTerm('GO:0000001', 'positive regulation of test process',
                  'biological_process', 'This is another test GO term.',
                  is_a=['GO:0000000'])
    return go_term


def test_basic(my_go_term, my_other_term):

    ontology = GeneOntology([my_go_term, my_other_term])

    assert isinstance(ontology, GeneOntology)
    assert isinstance(repr(ontology), str)
    assert isinstance(str(ontology), str)
    assert isinstance(text(ontology), text)
    assert isinstance(ontology.hash, text)

    # test access methods
    assert len(ontology) == 2
    assert my_go_term.id in ontology
    assert ontology[my_go_term.id] == my_go_term
    del ontology[my_go_term.id]
    assert my_go_term.id not in ontology
    ontology[my_go_term.id] = my_go_term
    assert my_go_term.id in ontology

    # test additional access methods
    assert ontology.get_term_by_id(my_go_term.id) == my_go_term
    assert ontology.get_term_by_acc(my_go_term.acc) == my_go_term

    # test comparisons
    other = copy.deepcopy(ontology)
    assert other == ontology
    del other[my_other_term.id]
    assert other != ontology

    # test iteration
    assert set(list(iter(ontology))) == set([my_go_term, my_other_term])


@pytest.mark.online
def test_real(my_gene_ontology):
    assert isinstance(my_gene_ontology, GeneOntology)