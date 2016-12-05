# Copyright (c) 2016 Florian Wagner
#
# This file is part of GenomeTools.
#
# GenomeTools is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License, Version 3,
# as published by the Free Software Foundation.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

"""Tests for the `GOAnnotation` class and the GO annotation file parser."""

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import str as text

import copy

import pytest

from genometools.ontology import GeneOntology, GOAnnotation, \
                                 parse_gaf, get_goa_gene_sets
from genometools.basic import GeneSetCollection

@pytest.fixture
def my_go_annotation(my_go_term):
    go_annotation = GOAnnotation(
        db='UniProtKB',
        db_id='A0A024R161',
        db_symbol='DNAJC25-GNG10',
        go_term=my_go_term,
        db_ref='GO_REF:0000002',
        ev_code='IEA',
        db_type='protein',
        taxon='taxon:9606',
        date='20161001',
        assigned_by='UniProt',
        
        qualifier=None,
        with_from='InterPro:IPR001770|InterPro:IPR015898',
        db_name='Guanine nucleotide-binding protein subunit gamma',
        db_syn=['A0A024R161_HUMAN', 'DNAJC25-GNG10', 'hCG_1994888']
    )
    return go_annotation


def test_basic(my_go_annotation):

    assert isinstance(my_go_annotation, GOAnnotation)
    assert isinstance(repr(my_go_annotation), str)
    assert isinstance(str(my_go_annotation), str)
    assert isinstance(text(my_go_annotation), text)
    assert isinstance(my_go_annotation.hash, text)
    
    # test comparisons
    other = copy.deepcopy(my_go_annotation)
    assert other == my_go_annotation
    other.ev_code = 'IDA'
    assert other != my_go_annotation


def test_list(my_go_annotation, my_go_term):

    gene_ontology = GeneOntology([my_go_term])
    l = my_go_annotation.to_list()
    assert isinstance(l, list)

    other = GOAnnotation.from_list(gene_ontology, l)
    assert isinstance(other, GOAnnotation)
    assert other == my_go_annotation


@pytest.mark.online
def test_parser(my_goa_file, my_genome, my_gene_ontology):

    ev_codes = ['IEA', 'TAS']
    go_annotations = parse_gaf(my_goa_file, my_gene_ontology,
                               genome=my_genome,
                               db='UniProtKB', ev_codes=ev_codes)
    
    assert isinstance(go_annotations, list)
    assert len(go_annotations) > 0
    for ann in go_annotations:
        assert isinstance(ann, GOAnnotation)
        assert ann.ev_code in ev_codes

    gene_sets = get_goa_gene_sets(go_annotations)
    assert isinstance(gene_sets, GeneSetCollection)