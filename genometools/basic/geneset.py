# Copyright (c) 2015, 2016 Florian Wagner
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

"""Module containing the `GeneSet` class.

"""

from collections import Iterable

class GeneSet(object):
    """A gene set."""

    def __init__(self, id_, name, genes, source = None,
            collection = None, description = None):

        assert isinstance(name, (str, unicode))
        assert isinstance(genes, Iterable)
        for g in genes:
            assert isinstance(g, (str, unicode))

        if id_ is not None:
            assert isinstance(id_, (str, unicode))
        if source is not None:
            assert isinstance(source, (str, unicode))
        if collection is not None:
            assert isinstance(collection, (str, unicode))
        if description is not None:
            assert isinstance(description, (str, unicode))

        self.id = id_
        self.name = name
        self.genes = frozenset(genes)
        self.source = source
        self.collection = collection
        self.description = description

    def __eq__(self):
        return '<%s (hash = %d)>' %(self.__class__.__name__, hash(self))

    def __str__(self):
        return '<%s with %d genes>' %(self.__class__.__name__, self.size)

    def __hash__(self):
        return hash(tuple(self.to_list()))

    @property
    def size(self):
        return len(self.genes)

    def to_list(self):
        l = []
        src = self.source or ''
        coll = self.collection or ''
        id_ = self.id or ''
        desc = self.description or ''

        l.append(id_)
        l.append(src)
        l.append(coll)
        l.append(self.name)
        l.append(','.join(sorted(self.genes)))
        l.append(desc)
        return l

    @classmethod
    def from_list(cls, l):

        id_ = l[0]
        name = l[3]
        genes = l[4].split(',')

        src = l[1] or None
        coll = l[2] or None
        desc = l[5] or None

        return cls(id_, name, genes, src, coll, desc)
