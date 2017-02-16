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

"""Tests for the `ExpGene` class."""

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import str as text

from genometools.expression import ExpGene


def test_init(my_genes):
    for g in my_genes:
        assert isinstance(g, ExpGene)
        assert isinstance(repr(g), str)
        assert isinstance(str(g), str)
        assert isinstance(text(g), text)


def test_list(my_genes):
    for g in my_genes:
        other = ExpGene.from_dict(g.to_dict())
        assert other is not g
        assert other == g