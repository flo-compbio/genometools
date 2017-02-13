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

"""Module containing the `ExpGene` class.

"""

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
_oldstr = str
from builtins import *

import logging

logger = logging.getLogger(__name__)


class ExpGene(object):
    """A gene in a gene expression analysis.

    Instances are to be treated as immutable, to allow use of ExpGene
    objects to be used in sets etc.

    Parameters
    ----------
    name : str
        See :attr:`name` attribute.
    chromosome : str or None, optional
        See :attr:`chromosome` attribute. [None]
    position : int or None, optional
        See :attr:`position` attribute. [None]
    length : int or None, optional
        See :attr:`length` attribute. [None] 
    ensembl_id : str or None, optional
        See :attr:`ensembl_id` attribute. [None]

    Attributes
    ----------
    name : str
        The gene name (use the official gene symbol, if available).
    chromosome : str or None
        The chromosome that the gene is located on.
    position : int or None
        The chromosomal location (base-pair index) of the gene.
        The sign of the this attribute indicates whether the gene is on the
        plus or minus strand. Base pair indices are 0-based.
    ensembl_id : list of str
        The Ensembl ID of the gene.
    """
    def __init__(self, name,
                 chromosome=None, position=None, length=None,
                 ensembl_id=None, source=None, type_=None):

        # type checks
        assert isinstance(name, (str, _oldstr))
        if chromosome is not None:
            assert isinstance(chromosome, (str, _oldstr))
        if position is not None:
            assert isinstance(position, int)
        if length is not None:
            assert isinstance(length, int)
        if ensembl_id is not None:
            assert isinstance(ensembl_id, (str, _oldstr))
        if source is not None:
            assert isinstance(source, (str, _oldstr))
        if type_ is not None:
            assert isinstance(type_, (str, _oldstr))

        self._name = name
        self._chromosome = chromosome
        self._position = position
        self._length = length
        self._ensembl_id = ensembl_id
        self.source = source
        self.type_ = type_

    def __repr__(self):
        return '<%s "%s">' \
               % (self.__class__.__name__, self._name)

    def __str__(self):
        return '<%s "%s" (Chromosome: %s, Position: %s, Length: %s, ' \
               'Ensembl ID: %s, Source: %s, Type: %s)>' \
               % (self.__class__.__name__, self._name,
                  self._chromosome, str(self._position), str(self._length),
                  self._ensembl_id, self.source, self.type_)

    def __eq__(self, other):
        if self is other:
            return True
        elif type(self) is type(other):
            return repr(self) == repr(other)
        else:
            return NotImplemented

    def __ne__(self, other):
        return not self.__eq__(other)

    def __lt__(self, other):
        if self is other:
            return False
        if type(self) is type(other):
            if self._name < other._name:
                return True
            else:
                return False
        else:
            return NotImplemented

    def __hash__(self):
        return hash(repr(self))

    @property
    def name(self):
        return self._name

    @property
    def chromosome(self):
        return self._chromosome
    
    @property
    def position(self):
        return self._position

    @property
    def length(self):
        return self._length

    @property
    def ensembl_id(self):
        return self._ensembl_id

    #def to_list(self):
    #    return [self._name, ','.join(self._chromosomes),
    #            ','.join(self._ensembl_ids)]

    def to_dict(self):
        return {
            'name': self.name,
            'chromosome': self.chromosome or '',
            'position': self.position or '',
            'length': self.length or '',
            'ensembl_id': self.ensembl_id or '',
            'source': self.source or '',
            'type': self.type_ or '',
        }

    @classmethod
    def from_dict(cls, data):
        """Generate an `ExpGene` gene object from a dictionary.

        Parameters
        ----------
        data : dict
            A dictionary with keys corresponding to attribute names.
            Attributes with missing keys will be assigned `None`.
            See also :meth:`to_list`.

        Returns
        -------
        `ExpGene`
            The gene.
        """
        assert isinstance(data, dict)
        assert 'name' in data  # required

        # make a copy
        data = dict(data)

        for attr in ['chromosome', 'ensembl_id', 'position', 'length',
                     'source', 'type']:
            if attr in data and data[attr] == '':
                data[attr] = None            

        data['type_'] = data['type']
        del data['type']

        return cls(**data)
