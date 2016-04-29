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

"""Module containing the `GeneSet` class."""

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import *


class GeneSet(object):
    """A gene set.

    A gene set is just what the name implies: A set of genes. Usually, gene
    sets are used to group genes that share a certain property (e.g., genes
    that perform related functions, or genes that are frequently co-expressed).
    The genes in the gene set are not ordered.
    
    Parameters
    ----------
    id: str
        See :attr:`id` attribute.
    name: str
        See :attr:`name` attribute.
    genes: set, list or tuple of str
        See :attr:`genes` attribute.
    source: str, optional
        See :attr:`source` attribute. (None)
    collection: str, optional
        See :attr:`collection` attribute. (None)
    description: str, optional
        See :attr:`description` attribute. (None)

    Attributes
    ----------
    id_: str
        The (unique) ID of the gene set.
    name: str
        The name of the gene set.
    genes: set of str
        The list of genes in the gene set.
    source: None or str
        The source / origin of the gene set (e.g., "MSigDB")
    collection: None or str
        The collection that the gene set belongs to (e.g., "c4" for gene sets
        from MSigDB).
    description: None or str
        The description of the gene set.
    """
    def __init__(self, id, name, genes,
                 source=None, collection=None, description=None):

        assert isinstance(id, str)
        assert isinstance(name, str)
        assert isinstance(genes, (set, list, tuple))
        for g in genes:
            assert isinstance(g, str)

        if source is not None:
            assert isinstance(source, str)
        if collection is not None:
            assert isinstance(collection, str)
        if description is not None:
            assert isinstance(description, str)

        self.id = id
        self.name = name
        self.genes = set(genes)
        self.source = source
        self.collection = collection
        self.description = description

    @property
    def _gene_str(self):
        return ', '.join('"%s"' % g for g in sorted(self.genes))

    @property
    def _source_str(self):
        return '"%s"' % self.source \
            if self.source is not None else 'None'

    @property
    def _coll_str(self):
        return '"%s"' % self.collection \
            if self.collection is not None else 'None'

    @property
    def _desc_str(self):
        return '"%s"' % self.description \
            if self.description is not None else 'None'

    def __repr__(self):
        return '<%s(id="%s", name="%s", genes=[%s], source=%s, ' \
               'collection=%s, description=%s)' \
                % (self.__class__.__name__, self.id, self.name, self._gene_str,
                   self._source_str, self._coll_str, self._desc_str)

    def __str__(self):
        return '<%s "%s" (id=%s; source=%s; collection=%s; size=%d)>' \
                % (self.__class__.__name__, self.name,
                   self.id, self._source_str, self._coll_str, self.size)

    def __eq__(self, other):
        if self is other:
            return True
        elif type(self) is type(other):
            return self.__dict__ == other.__dict__
        else:
            return NotImplemented

    def __ne__(self, other):
        return not self.__eq__(other)

    @property
    def size(self):
        """The size of the gene set (i.e., the number of genes in it)."""
        return len(self.genes)

    def to_list(self):
        """Converts the GeneSet object to a flat list of strings.

        Note: see also :meth:`from_list`.

        Parameters
        ----------

        Returns
        -------
        list of str
            The data from the GeneSet object as a flat list.
        """
        src = self.source or ''
        coll = self.collection or ''
        desc = self.description or ''

        l = [self.id, src, coll, self.name, ','.join(sorted(self.genes)), desc]
        return l

    @classmethod
    def from_list(cls, l):
        """Generate an GeneSet object from a list of strings.

        Note: See also :meth:`to_list`.

        Parameters
        ----------
        l: list or tuple of str
            A list of strings representing gene set ID, name, genes,
            source, collection, and description. The genes must be
            comma-separated. See also :meth:`to_list`.

        Returns
        -------
        `genometools.basic.GeneSet`
            The gene set.
        """
        id_ = l[0]
        name = l[3]
        genes = l[4].split(',')

        src = l[1] or None
        coll = l[2] or None
        desc = l[5] or None

        return cls(id_, name, genes, src, coll, desc)
