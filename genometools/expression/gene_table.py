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

"""Module containing the `ExpGeneTable` class."""

import logging
import hashlib
from typing import List, Dict

import pandas as pd
from pandas.util import hash_pandas_object

from .. import ensembl
from .gene import ExpGene

_LOGGER = logging.getLogger(__name__)


class ExpGeneTable(pd.DataFrame):
    """A table with the complete set of genes in a gene expression analysis.

    Parameters
    ----------
    See `pandas.DataFrame` constructor.

    Attributes
    ----------
    See `pandas.DataFrame` attributes.
    
    Notes
    -----
    This class inherits from `pandas.DataFrame`, and just adds
    some basic functionalities to it.
    """
    def __init__(self, *args, **kwargs):
        # call base class constructor
        pd.DataFrame.__init__(self, *args, **kwargs)

        if self.index.name is None:
            self.index.name = 'ensembl_id'

        # make sure columns are ordered correctly
        cols = ['name', 'chromosome', 'position', 'length', 'type', 'source']
        if self.columns.tolist() != cols:
            self = self[cols]

    @property
    def _constructor(self):
        return ExpGeneTable
    
    def __repr__(self):
        return '<%s object (%d genes, hash="%s")>' \
                % (self.__class__.__name__, len(self), self.hash)

    def __str__(self):
        types = ', '.join('%s (%d)' %(t, c)
                          for t, c in self['type'].value_counts().iteritems())
        return '<%s with %d genes -- types: %s>\n%s' \
               % (self.__class__.__name__, len(self),
                  types, pd.DataFrame.__str__(self))

    def __contains__(self, gene_or_id):
        """Tests whether a gene or ID is present."""
        if isinstance(gene_or_id, str):
            # test if Ensembl ID is contained
            return gene_or_id in self.index
        elif isinstance(gene_or_id, ExpGene):
            return gene_or_id.ensembl_id in self.index
        else:
            return NotImplemented

    @property
    def hash(self):
        """Generate a hash value."""
        h = hash_pandas_object(self, index=True)
        return hashlib.md5(h.values.tobytes()).hexdigest()
    
    @property
    def num_genes(self):
        """Return the number of genes."""
        return len(self.index)

    @property
    def genes(self):
        """Return a list of all genes."""
        return [ExpGene.from_series(g)
                for i, g in self.reset_index().iterrows()]
    
    @property
    def gene_names(self):
        """Return a list of all gene names."""
        return self['name'].tolist()

    @property
    def gene_ids(self):
        """Return a list of all gene IDs."""
        return self.index.tolist()
    
    @classmethod
    def read_tsv(cls, file_or_buffer: str):
        """Read genes from tab-delimited text file."""
        df = pd.read_csv(file_or_buffer, sep='\t', index_col=0)
        df = df.where(pd.notnull(df), None)
        # Note: df.where(..., None) changes all column types to `object`.
        return cls(df)
    
    @classmethod
    def from_genes(cls, genes: List[ExpGene]):
        """Initialize instance using a list of `ExpGene` objects."""
        data = [g.to_dict() for g in genes]
        index = [d.pop('ensembl_id') for d in data]
        table = cls(data, index=index)
        return table
    
    @classmethod
    def from_gene_ids(cls, gene_ids: List[str]):
        """Initialize instance from gene IDs."""
        genes = [ExpGene(id_) for id_ in gene_ids]
        return cls.from_genes(genes)

    @classmethod
    def from_gene_ids_and_names(cls, gene_names: Dict[str, str]):
        """Initialize instance from gene IDs and names."""
        genes = [ExpGene(id_, name=name) for id_, name in gene_names.items()]
        return cls.from_genes(genes)

    @classmethod
    def read_gtf(cls, file_or_buffer, **kwargs):
        genes = ensembl.get_protein_coding_genes(file_or_buffer, **kwargs)
        return cls(genes)

    def write_tsv(self, file_or_buffer: str):
        """Write genes to tab-delimited text file."""
        self.to_csv(file_or_buffer, sep='\t')

    def find_gene(self, name: str):
        """Find gene(s) by name."""
        result = [ExpGene.from_series(s)
                  for i, s in self.loc[self['name'] == name].iterrows()]
        return result
