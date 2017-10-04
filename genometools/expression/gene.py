# Copyright (c) 2015-2017 Florian Wagner
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

import logging
from typing import Dict, Union

import pandas as pd

_LOGGER = logging.getLogger(__name__)


class ExpGene:
    """A gene in a gene expression analysis.

    Instances are to be treated as immutable, to allow ExpGene
    objects to be used in sets etc.

    Parameters
    ----------
    ensembl_id : str
        See :attr:`ensembl_id` attribute.
    name : str, optional
        See :attr:`name` attribute. [None]
    chromosome : str or None, optional
        See :attr:`chromosome` attribute. [None]
    position : int or None, optional
        See :attr:`position` attribute. [None]
    length : int or None, optional
        See :attr:`length` attribute. [None]
    type_ : str, optional
        See :attr:`type` attribute. [None]
    source : str, optional
        See :attr:`source` attribute. [None]

    Attributes
    ----------
    ensembl_id : str
        The Ensembl ID of the gene.
    name : str or None
        The gene name (sometimes called "gene symbol").
    chromosome : str or None
        The chromosome that the gene is located on.
    position : int or None
        The chromosomal location of the gene, defined as the index of
        its 5'-most base, according to its orientation on the chromosome.
        The sign of the value indicates whether the gene is located on the
        plus or minus strand.
    length : int or None
        The length of the gene.
    type : str or None
        The type ("biotype") of the gene, using Ensembl identifiers.
        Examples: "protein_coding", "polymorphic_pseudogene".
    source : str or None
        The source of the gene, using Ensembl identifiers.
        Examples: "havana", "ensembl_havana".
    """
    def __init__(self, ensembl_id: str, name: str = None,
                 chromosome: str = None, position: int = None,
                 length: int = None,
                 type_: str = None, source: str = None):

        self._ensembl_id = ensembl_id
        self._name = name
        self._chromosome = chromosome
        self._position = position
        self._length = length
        self._ensembl_id = ensembl_id
        self._type = type_
        self._source = source

    def __repr__(self):
        # there should only be one ExpGene object per Ensembl ID
        return '<%s %s>' \
               % (self.__class__.__name__, self._ensembl_id)

    def __str__(self):
        name = self._name or '[no name]'
        type_ = self._type or '[unknown]'
        return '<%s %s (%s) of type "%s" ' \
               '(Chromosome: %s, Position: %s, Length: %s, Source: %s)>' \
               % (self.__class__.__name__, self._ensembl_id, name, type_,
                  self._chromosome, str(self._position), str(self._length),
                  self._source)

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
            if self._ensembl_id < other._ensembl_id:
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
    
    @property
    def source(self):
        return self._source
    
    @property
    def type(self):
        return self._type

    @classmethod
    def from_dict(cls, data: Dict[str, Union[str, int]]):
        """Generate an `ExpGene` object from a dictionary.

        Parameters
        ----------
        data : dict
            A dictionary with keys corresponding to attribute names.
            Attributes with missing keys will be assigned `None`.

        Returns
        -------
        `ExpGene`
            The gene.
        """
        assert isinstance(data, dict)

        if 'ensembl_id' not in data:
            raise ValueError('An "ensembl_id" key is missing!')

        # make a copy
        data = dict(data)
        
        for attr in ['name', 'chromosome', 'position', 'length',
                     'type', 'source']:
            if attr in data and data[attr] == '':
                data[attr] = None            

        data['type_'] = data['type']
        del data['type']

        return cls(**data)

    @classmethod
    def from_series(cls, s: pd.Series):
        """Generate an `ExpGene` object from a pandas Series."""
        return cls.from_dict(s.to_dict())

    def to_dict(self):
        d = dict([attr, getattr(self, attr)] for attr in
                 ['ensembl_id', 'name', 'chromosome', 'position', 'length',
                  'type', 'source'])
        return d
