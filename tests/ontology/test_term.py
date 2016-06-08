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

from genometools.ontology import GOTerm

@pytest.fixture
def my_term():
    term = GOTerm('GO:0000000', 'positive regulation of test process',
                  'biological_process', 'This is a test GO term.')
    return term

def test_basic(my_term):
    assert isinstance(my_term, GOTerm)
    assert isinstance(repr(my_term), str)
    assert isinstance(str(my_term), str)
    assert isinstance(text(my_term), text)

    assert isinstance(my_term.id, text)
    assert isinstance(my_term.name, text)
    assert isinstance(my_term.domain, text)
    assert isinstance(my_term.definition, text)

    assert isinstance(my_term.is_a, set)
    assert isinstance(my_term.part_of, set)
    assert isinstance(hash(my_term), int)

def test_compare(my_term):
    """Test comparison between two GOTerm objects."""
    other = copy.deepcopy(my_term)
    assert other == my_term

    other.name = 'new name'
    assert other == my_term

    other.id = 'GO:0000001'
    assert other != my_term

def test_pretty_format(my_term):
    # call get_pretty_format with all arguments
    fmt = my_term.get_pretty_format(include_id=True, max_name_length=8,
                                    abbreviate=True)
    assert isinstance(fmt, text)  # make sure return type is correct
    assert 'pos.' in fmt  # make sure abbreviations work
    # make sure name shortening works
    assert len(fmt) == len('BP: [] (GO:1234567)') - 2 + 8
