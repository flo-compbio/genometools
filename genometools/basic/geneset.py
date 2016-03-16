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

class GeneSet(object):
    """A gene set.

    A gene set is just what the name implies: A set of genes. Usually, gene
    sets are used to group genes that share a certain property (e.g., genes
    that perform related functions, or genes that are frequently co-expressed).
    The genes in the gene set are not ordered.
    
    Parameters
    ----------
    id_: str or unicode
        See :attr:`id` attribute.
    name: str or unicode
        See :attr:`name` attribute.
    genes: set, frozenset, list or tuple of (str or unicode)
        See :attr:`genes` attribute.
    source: str or unicode, optional
        See :attr:`source` attribute. (None)
    collection: str or unicode, optional
        See :attr:`collection` attribute. (None)
    description: str or unicode, optional
        See :attr:`description` attribute. (None)

    Attributes
    ----------
    id_: str or unicode
        The (unique) ID of the gene set.
    name: str or unicode
        The name of the gene set.
    genes: frozenset of (str or unicode)
        The list of genes in the gene set.
    source: None or (str or unicode)
        The source / origin of the gene set (e.g., "MSigDB")
    collection: None or (str or unicode)
        The collection that the gene set belongs to (e.g., "c4" for gene sets
        from MSigDB).
    description: None or (str or unicode)
        The description of the gene set.
    """

    def __init__(self, id_, name, genes, source = None,
            collection = None, description = None):

        assert isinstance(id_, (str, unicode))
        assert isinstance(name, (str, unicode))
        assert isinstance(genes, (list, tuple, set, frozenset))
        for g in genes:
            assert isinstance(g, (str, unicode))

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

    def __repr__(self):

        src_str = unicode(self.source)
        coll_str = unicode(self.collection)
        desc_str = unicode(self.description)

        return '<%s "%s" (id=%s; source=%s; collection=%s; size=%d; hash=%d)>' \
                %(self.__class__.__name__, self.name,
                    self.id, self.source, self.collection,
                    self.size, hash(self))

    def __str__(self):
        return '<%s "%s" (id=%s; source=%s; collection=%s; size=%d)>' \
                %(self.__class__.__name__, self.name,
                    self.id, self.source, self.collection, self.size)

    def __eq__(self, other):
        if self is other:
            return True
        elif type(self) != type(other):
            return False
        else:
            return repr(self) == repr(other)

    def __ne__(self, other):
        return not (self == other)

    def __hash__(self):
        data = []
        data.append(self.id)
        data.append(self.name)
        data.append(self.genes)
        data.append(self.source)
        data.append(self.collection)
        data.append(self.description)
        return hash(tuple(data))

    @property
    def size(self):
        """The size of the gene set (i.e., the number of genes in it)."""
        return len(self.genes)

    def to_list(self):
        """Converts the GeneSet object to a flat list of strings.

        Note: see also :meth:`from_list`.

        Parameters
        ----------
        None
        
        Returns
        -------
        list of str
            The data from the GeneSet object as a flat list.
        """
        l = []
        src = self.source or ''
        coll = self.collection or ''
        desc = self.description or ''

        l.append(self.id)
        l.append(src)
        l.append(coll)
        l.append(self.name)
        l.append(','.join(sorted(self.genes)))
        l.append(desc)
        return l

    @classmethod
    def from_list(cls, l):
        """Generate an GeneSet object from a list of strings.

        Note: See also :meth:`to_list`.

        Parameters
        ----------
        l: list or tuple of (str or unicode)
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
