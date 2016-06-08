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

"""Module containing the `GeneOntology` class.
"""

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import *

import gzip
import hashlib
# import re
import six
# import sys
import logging
# import bisect

from collections import OrderedDict, Iterable

import unicodecsv as csv

from genometools import misc
from . import GOTerm

if six.PY2:
    import cPickle as pickle
else:
    import pickle

logger = logging.getLogger(__name__)


class GeneOntology(object):
    """A Gene Ontology.

    This class provides functions for parsing text files describing the Gene
    Ontology, and for accessing information about specific GO terms.

    Parameters
    ----------

    Attributes
    ----------

    Methods
    -------
    get_term_by_id(id_)
        Return the term with the given term ID as a `GOTerm` object.
    get_term_by_name(name)
        Return the term with the given name as a `GOTerm` object.
    save(ofn, compress=False)
        Stores the GOParser object as a `pickle` file. If ``compress`` is set
        to True, the object is stored as a gzip'ed pickle file.
    load(fn)
        Loads the `GOParser` object from a `pickle` file. Gzip compression is
        detected automatically.

    Examples
    --------
    The following example assumes that the Gene Ontology OBO file has been downloaded.
    >>> from genometools.annotation import GeneOntology
    >>> ontology = GeneOntology.read_obo('go-basic.obo')
    """
    # TODO: finish docstring
    def __init__(self, terms=None, syn2id=None, alt_id=None, name2id=None):

        if terms is None:
            terms = []

        if syn2id is None:
            syn2id = {}

        if alt_id is None:
            alt_id = {}

        if name2id is None:
            name2id = {}

        assert isinstance(terms, Iterable)
        assert isinstance(syn2id, dict)
        assert isinstance(alt_id, dict)
        assert isinstance(name2id, dict)

        term_dict = {}
        for t in terms:
            term_dict[t.id] = t
        self._term_dict = term_dict

        self.syn2id = syn2id
        self.alt_id = alt_id
        self.name2id = name2id
        self._flattened = False

    def __repr__(self):
        return '<%s instance (%d GO terms, hash="%s")>' \
               % (self.__class__.__name__, len(self), self.hash)

    def __str__(self):
        return '<%s instance with %d GO terms>' \
               % (self.__class__.__name__, len(self))

    def __eq__(self, other):
        if self is other:
            return True
        elif type(self) is type(other):
            return repr(self) == repr(other)
        else:
            raise NotImplemented

    def __ne__(self, other):
        return not self.__eq__(other)

    def __getitem__(self, key):
        return self._term_dict[key]

    def __setitem__(self, key, value):
        assert isinstance(value, GOTerm)
        self._term_dict[key] = value

    def __delitem__(self, key):
        del self._term_dict[key]

    def __len__(self):
        return len(self._term_dict)

    def __contains__(self, key):
        return key in self._term_dict

    def __iter__(self):
        return iter(self._term_dict.values())

    @property
    def hash(self):
        data_str = ';'.join(
            [repr(self._term_dict[id_]) for id_ in sorted(self._term_dict.keys())] +
            [repr(var) for var in [self.syn2id, self.alt_id, self.name2id]]
        )
        data = data_str.encode('UTF-8')
        return str(hashlib.md5(data).hexdigest())

    @property
    def flattened(self):
        return self._flattened

    def write_pickle(self, path, compress=False):
        """Serialize the current `GOParser` object and store it in a pickle file.

        Parameters
        ----------
        path: str
            Path of the output file.
        compress: bool, optional
            Whether to compress the file using gzip.

        Returns
        -------
        None

        Notes
        -----
        Compression with gzip is significantly slower than storing the file
        in uncompressed form.
        """
        logger.info('Writing pickle to "%s"...', path)
        if compress:
            with gzip.open(path, 'wb') as ofh:
                pickle.dump(self, ofh, pickle.HIGHEST_PROTOCOL)
        else:
            with open(path, 'wb') as ofh:
                pickle.dump(self, ofh, pickle.HIGHEST_PROTOCOL)

    @staticmethod
    def read_pickle(fn):
        """Load a GOParser object from a pickle file.

        The function automatically detects whether the file is compressed
        with gzip.

        Parameters
        ----------
        fn: str
            Path of the pickle file.

        Returns
        -------
        `GOParser`
            The GOParser object stored in the pickle file.
        """
        with misc.open_plain_or_gzip(fn, 'rb') as fh:
            parser = pickle.load(fh)
        return parser

    def get_term_by_id(self, id_):
        """Get the GO term corresponding to the given GO term ID.

        Parameters
        ----------
        id_: str
            A GO term ID.

        Returns
        -------
        `GOTerm`
            The GO term corresponding to the given ID.
        """
        return self[id_]

    def get_term_by_acc(self, acc):
        """Get the GO term corresponding to the given GO term accession number.

        Parameters
        ----------
        acc: int
            The GO term accession number.

        Returns
        -------
        `GOTerm`
            The GO term corresponding to the given accession number.
        """
        return self[GOTerm.acc2id(acc)]

    def get_term_by_name(self, name):
        """Get the GO term with the given GO term name.

        If the given name is not associated with any GO term, the function will
        search for it among synonyms.

        Parameters
        ----------
        name: str
            The name of the GO term.

        Returns
        -------
        `GOTerm`
            The GO term with the given name.

        Raises
        ------
        ValueError
            If the given name is found neither among the GO term names, nor
            among synonyms.
        """
        term = None
        try:
            term = self.terms[self.name2id[name]]
        except KeyError:
            try:
                term = self.terms[self.syn2id[name]]
            except KeyError:
                pass
            else:
                logger.info('GO term name "%s" is a synonym for "%s".',
                            name, term.name)

        if term is None:
            raise ValueError('GO term name "%s" not found!' % name)

        return term

    @classmethod
    def read_obo(cls, path, flatten=True, part_of_cc_only=False):
        """ Parse an OBO file and store GO term information.

        Parameters
        ----------
        path: str
            Path of the OBO file.
        flatten: bool, optional
            If set to False, do not generate a list of all ancestors and
            descendants for each GO term.
        part_of_cc_only: bool, optional
            Legacy parameter for backwards compatibility. If set to True,
            ignore ``part_of`` relations outside the ``cellular_component``
            domain.

        Notes
        -----
        The OBO file must end with a line break.
        """

        name2id = {}
        alt_id = {}
        syn2id = {}
        terms = []

        with open(path) as fh:
            n = 0
            while True:
                try:
                    nextline = next(fh)
                except StopIteration:
                    break
                if nextline == '[Term]\n':
                    n += 1
                    id_ = next(fh)[4:-1]
                    # acc = get_acc(id_)
                    name = next(fh)[6:-1]
                    name2id[name] = id_
                    domain = next(fh)[11:-1]
                    def_ = None
                    is_a = set()
                    part_of = set()
                    l = next(fh)
                    while l != '\n':
                        if l.startswith('alt_id:'):
                            alt_id[l[8:-1]] = id_
                        elif l.startswith('def: '):
                            idx = l[6:].index('"')
                            def_ = l[6:(idx+6)]
                        elif l.startswith('is_a:'):
                            is_a.add(l[6:16])
                        elif l.startswith('synonym:'):
                            idx = l[10:].index('"')
                            if l[(10+idx+2):].startswith("EXACT"):
                                s = l[10:(10+idx)]
                                syn2id[s] = id_
                        elif l.startswith('relationship: part_of'):
                            if part_of_cc_only:
                                if domain == 'cellular_component':
                                    part_of.add(l[22:32])
                            else:
                                part_of.add(l[22:32])
                        l = next(fh)
                    assert def_ is not None
                    terms.append(GOTerm(id_, name, domain, def_, is_a, part_of))

        logger.info('Parsed %d GO term definitions.', n)

        ontology = cls(terms, syn2id, alt_id, name2id)

        # store children and parts
        logger.info('Adding child and part relationships...')
        for term in ontology:
            for parent in term.is_a:
                ontology[parent].children.add(term.id)
            for whole in term.part_of:
                ontology[whole].parts.add(term.id)

        if flatten:
            logger.info('Flattening ancestors...')
            ontology._flatten_ancestors()
            logger.info('Flattening descendants...')
            ontology._flatten_descendants()
            ontology._flattened = True

        return ontology

    def _flatten_ancestors(self, include_part_of=True):
        """Determines and stores all ancestors of each GO term.

        Parameters
        ----------
        include_part_of: bool, optional
            Whether to include ``part_of`` relations in determining
            ancestors.

        Returns
        -------
        None
        """
        def get_all_ancestors(term):
            ancestors = set()
            for id_ in term.is_a:
                ancestors.add(id_)
                ancestors.update(get_all_ancestors(self[id_]))
            if include_part_of:
                for id_ in term.part_of:
                    ancestors.add(id_)
                    ancestors.update(get_all_ancestors(self[id_]))
            return ancestors

        for term in self:
            term.ancestors = get_all_ancestors(term)

    def _flatten_descendants(self, include_parts=True):
        """Determines and stores all descendants of each GO term.

        Parameters
        ----------
        include_parts: bool, optional
            Whether to include ``part_of`` relations in determining
            descendants.

        Returns
        -------
        None
        """
        def get_all_descendants(term):
            descendants = set()
            for id_ in term.children:
                descendants.add(id_)
                descendants.update(get_all_descendants(self[id_]))
            if include_parts:
                for id_ in term.parts:
                    descendants.add(id_)
                    descendants.update(get_all_descendants(self[id_]))
            return descendants

        for term in self:
            term.descendants = get_all_descendants(term)