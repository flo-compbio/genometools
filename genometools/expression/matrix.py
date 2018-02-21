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

"""Module containing the `ExpMatrix` class."""

import logging
import importlib
import hashlib
from collections import OrderedDict
from typing import Iterable
import tarfile
import io

import pandas as pd
import numpy as np
import scipy as sp
import scipy.io
from scipy import sparse
import csv
#import six

# from .. import misc
from . import ExpGene, ExpGeneTable
profile = importlib.import_module('.profile', package='genometools.expression')
# - "import profile" is not possible, since a "profile" module exists
#    in the standard library
# - "from . import profile" fails due to cyclical imports

_LOGGER = logging.getLogger(__name__)


class ExpMatrix(pd.DataFrame):
    """A gene expression matrix.

    This class inherits from `pandas.DataFrame`.

    Parameters
    ----------
    X : 2-dimensional `numpy.ndarray`
        See :attr:`X` attribute.

    Keyword-only Parameters
    -----------------------
    genes : list or tuple of str
        See :attr:`genes` attribute.
    samples : list or tuple of str
        See :attr:`samples` attribute.

    Additional Parameters
    ---------------------
    All `pandas.DataFrame` parameters.

    Attributes
    ----------
    genes : tuple of str
        The names of the genes (rows) in the matrix.
    samples : tuple of str
        The names of the samples (columns) in the matrix.
    X : 2-dimensional `numpy.ndarray`
        The matrix of expression values.
    """
    def __init__(self, *args,
                 X: np.ndarray = None,
                 genes: Iterable[str] = None, gene_label: str = None,
                 cells: Iterable[str] = None, cell_label: str = None,
                 samples: Iterable[str] = None, sample_label: str = None,
                 **kwargs):
        
        # check if user provided "X" keyword argument...
        # if so, pass it as the "data" argument to the DataFrame constructor
        if X is not None:
            if 'data' in kwargs:
                raise ValueError(
                    'Cannot specify both "X" and "data" arguments.')
            kwargs['data'] = X

        if cells is not None and samples is not None:
            raise ValueError('Cannot specify both "cells" and "samples" '
                             'arguments.')

        if cell_label is not None and sample_label is not None:
            raise ValueError('Cannot specify both "cell_label" and '
                             '"sample_label" arguments.')

        # call base class constructor
        pd.DataFrame.__init__(self, *args, **kwargs)

        # set gene (row) and cell/sample (column) labels, if provided
        if genes is not None:
            self.index = genes

        if cells is not None:
            self.columns = cells
        elif samples is not None:
            self.columns = samples

        # set index name (default: "Genes")
        if gene_label is not None:
            self.index.name = gene_label
        elif self.index.name is None:
            self.index.name = 'Genes'

        # set column name (default: "Cells")
        if cell_label is not None:
            self.columns.name = cell_label
        elif sample_label is not None:
            self.columns.name = sample_label
        elif self.columns.name is None:
            self.columns.name = 'Cells'


    def __eq__(self, other):
        if self is other:
            return True
        elif type(self) is type(other):
            return (self.index.equals(other.index) and \
                    self.columns.equals(other.columns) and \
                    self.equals(other))
        else:
            return pd.DataFrame.__eq__(self, other)

    def __ne__(self, other):
        return not self.__eq__(other)

    def __repr__(self):
        return '<%s instance (p=%d, n=%d, hash="%s">' \
               % (self.__class__.__name__, self.p, self.n, self.hash)

    #def __str__(self):
    #    return '<%s object with p=%d genes and n=%d samples>' \
    #           % (self.__class__.__name__, self.p, self.n)

    @property
    def hash(self):
        # warning: involves copying all the data
        gene_str = ','.join(str(s) for s in self.genes)
        sample_str = ','.join(str(s) for s in self.samples)
        data_str = ';'.join([gene_str, sample_str]) + ';'
        data = data_str.encode('UTF-8') + self.X.tobytes()
        return str(hashlib.md5(data).hexdigest())

    @property
    def _constructor(self):
        return ExpMatrix
    
    @property
    def _constructor_sliced(self):
        return profile.ExpProfile

    @property
    def p(self):
        """The number of genes."""
        return self.shape[0]

    @property
    def n(self):
        """The number of samples."""
        return self.shape[1]

    @property
    def genes(self):
        """Alias for `DataFrame.index`."""
        # return tuple(str(g) for g in self.index)
        return self.index

    @genes.setter
    def genes(self, gene_list):
        self.index = gene_list

    @property
    def cells(self):
        """Alias for `DataFrame.columns`."""
        # return tuple(str(s) for s in self.columns)
        return self.columns

    @cells.setter
    def cells(self, cell_list):
        self.columns = cell_list

    @property
    def samples(self):
        """Alias for `DataFrame.columns`."""
        # return tuple(str(s) for s in self.columns)
        return self.columns

    @samples.setter
    def samples(self, sample_list):
        self.columns = sample_list

    @property
    def X(self):
        """Alias for `DataFrame.values`."""
        return self.values

    def sum(self, axis=0, dtype=None):
        """Sum over rows or columns.
        
        Overrides the default `pandas.DataFrame.sum()` function to prevent
        temporary in-memory copies.
        """
        if axis not in [0, 1, 'index', 'columns']:
            raise ValueError('"axis" parameter must be one of 0, 1, "index", '
                             'or "columns".')

        sum_kwargs = {}
        if dtype is not None:
            sum_kwargs['dtype'] = dtype
        y = self.values.sum(axis=axis, **sum_kwargs)

        if axis == 0 or axis == 'index':
            y = profile.ExpProfile(y, genes=self.cells)
        else:
            y = profile.ExpProfile(y, genes=self.genes)

        return y


    def filter_variance(self, top):
        from .filter import filter_variance
        if top > self.p:
            raise ValueError('filter_variance() called with top=%d, but '
                             'there are only %d genes in the matrix.'
                             %(top, self.p))
        return filter_variance(self, top)
        

    def get_heatmap(self, highlight_genes=None, highlight_samples=None,
                          highlight_color=None, **kwargs):
        """Generate a heatmap (`ExpHeatmap`) of the matrix.

        See :class:`ExpHeatmap` constructor for keyword arguments.

        Parameters
        ----------
        highlight_genes : list of str
            List of genes to highlight
        highlight_color : str
            Color to use for highlighting

        Returns
        -------
        `ExpHeatmap`
            The heatmap.
        """
        from .visualize import ExpHeatmap
        from .visualize import HeatmapGeneAnnotation
        from .visualize import HeatmapSampleAnnotation

        #if highlight_genes is not None:
        #    assert isinstance(highlight_genes, Iterable)
        #if highlight_samples is not None:
        #    assert isinstance(highlight_genes, Iterable)
        if highlight_color is not None:
            assert isinstance(highlight_color, str)

        if highlight_color is None:
            highlight_color = 'blue'

        if highlight_genes is None:
            highlight_genes = []

        if highlight_samples is None:
            highlight_samples = []

        gene_annotations = kwargs.pop('gene_annotations', [])
        for g in highlight_genes:
            gene_annotations.append(
                HeatmapGeneAnnotation(g, highlight_color, label=g))

        sample_annotations = kwargs.pop('sample_annotations', [])
        for s in highlight_samples:
            sample_annotations.append(
                HeatmapSampleAnnotation(s, highlight_color, label=s)
            )

        return ExpHeatmap(self,
                          gene_annotations=gene_annotations,
                          sample_annotations=sample_annotations,
                          **kwargs)


    def get_figure(self, heatmap_kw=None, **kwargs):
        """Generate a plotly figure showing the matrix as a heatmap.

        This is a shortcut for ``ExpMatrix.get_heatmap(...).get_figure(...)``.

        See :func:`ExpHeatmap.get_figure` for keyword arguments.

        Parameters
        ----------
        heatmap_kw : dict or None
            If not None, dictionary containing keyword arguments to be passed
            to the `ExpHeatmap` constructor.

        Returns
        -------
        `plotly.graph_objs.Figure`
            The plotly figure.
        """
        if heatmap_kw is not None:
            assert isinstance(heatmap_kw, dict)

        if heatmap_kw is None:
            heatmap_kw = {}

        return self.get_heatmap(**heatmap_kw).get_figure(**kwargs)


    def sort_genes(self, stable=True, inplace=False, ascending=True):
        """Sort the rows of the matrix alphabetically by gene name.

        Parameters
        ----------
        stable: bool, optional
            Whether to use a stable sorting algorithm. [True]
        inplace: bool, optional
            Whether to perform the operation in place.[False]
        ascending: bool, optional
            Whether to sort in ascending order [True]
        
        Returns
        -------
        `ExpMatrix`
            The sorted matrix.
        """
        kind = 'quicksort'
        if stable:
            kind = 'mergesort'
        return self.sort_index(kind=kind, inplace=inplace, ascending=ascending)

    def sort_samples(self, stable=True, inplace=False, ascending=True):
        """Sort the columns of the matrix alphabetically by sample name.

        Parameters
        ----------
        stable: bool, optional
            Whether to use a stable sorting algorithm. [True]
        inplace: bool, optional
            Whether to perform the operation in place.[False]
        ascending: bool, optional
            Whether to sort in ascending order [True]

        Returns
        -------
        `ExpMatrix`
            The sorted matrix.
        """
        kind = 'quicksort'
        if stable:
            kind = 'mergesort'
        return self.sort_index(axis=1, kind=kind, inplace=inplace,
                               ascending=ascending)

    def center_genes(self, use_median=False, inplace=False):
        """Center the expression of each gene (row)."""
        if use_median:
            X = self.X - \
                np.tile(np.median(self.X, axis=1), (self.n, 1)).T
        else:
            X = self.X - \
                np.tile(np.mean(self.X, axis=1), (self.n, 1)).T

        if inplace:
            self.X[:,:] = X
            matrix = self
        else:
            matrix = ExpMatrix(genes=self.genes, samples=self.samples,
                               X=X)
        return matrix

    def standardize_genes(self, inplace=False):
        """Standardize the expression of each gene (row)."""
        matrix = self.center_genes(inplace=inplace)
        matrix.X[:,:] = matrix.X / \
            np.tile(np.std(matrix.X, axis=1, ddof=1), (matrix.n, 1)).T
        return matrix

    def filter_genes(self, gene_names : Iterable[str], inplace=False):
        """Filter the expression matrix against a _genome (set of genes).

        Parameters
        ----------
        gene_names: list of str
            The genome to filter the genes against.
        inplace: bool, optional
            Whether to perform the operation in-place.

        Returns
        -------
        ExpMatrix
            The filtered expression matrix.
        """

        return self.drop(set(self.genes) - set(gene_names),
                         inplace=inplace)

    @property
    def sample_correlations(self):
        """Returns an `ExpMatrix` containing all pairwise sample correlations.

        Returns
        -------
        `ExpMatrix`
            The sample correlation matrix.

        """
        C = np.corrcoef(self.X.T)
        corr_matrix = ExpMatrix(genes=self.samples, samples=self.samples, X=C)
        return corr_matrix


    @classmethod
    def read_tsv(cls, file_path: str, gene_table: ExpGeneTable = None,
                 encoding: str = 'UTF-8', sep: str = '\t'):
        """Read expression matrix from a tab-delimited text file.

        Parameters
        ----------
        file_path: str
            The path of the text file.
        gene_table: `ExpGeneTable` object, optional
            The set of valid genes. If given, the genes in the text file will
            be filtered against this set of genes. (None)
        encoding: str, optional
            The file encoding. ("UTF-8")
        sep: str, optional
            The separator. ("\t")

        Returns
        -------
        `ExpMatrix`
            The expression matrix.
        """
        # use pd.read_csv to parse the tsv file into a DataFrame
        matrix = cls(pd.read_csv(file_path, sep=sep, index_col=0, header=0,
                                 encoding=encoding))

        # parse index column separately
        # (this seems to be the only way we can prevent pandas from converting
        #  "nan" or "NaN" to floats in the index)['1_cell_306.120', '1_cell_086.024', '1_cell_168.103']
        #ind = pd.read_csv(file_path, sep=sep, usecols=[0, ], header=0,
        #                  encoding=encoding, na_filter=False)
        ind = pd.read_csv(file_path, sep=sep, usecols=[0, ], header=None,
                          skiprows=1, encoding=encoding, na_filter=False)

        matrix.index = ind.iloc[:, 0]
        matrix.index.name = 'Genes'

        if gene_table is not None:
            # filter genes
            matrix = matrix.filter_genes(gene_table.gene_names)

        return matrix


    def write_tsv(self, file_path: str, encoding: str = 'UTF-8',
                  sep: str = '\t'):
        """Write expression matrix to a tab-delimited text file.

        Parameters
        ----------
        file_path: str
            The path of the output file.
        encoding: str, optional
            The file encoding. ("UTF-8")

        Returns
        -------
        None
        """
        #if six.PY2:
        #    sep = sep.encode('UTF-8')

        self.to_csv(
            file_path, sep=sep, float_format='%.5f', mode='w',
            encoding=encoding, quoting=csv.QUOTE_NONE
        )

        _LOGGER.info('Wrote %d x %d expression matrix to "%s".',
                    self.p, self.n, file_path)


    def write_sparse(self, file_path: str):
        """Write a sparse representation to a tab-delimited text file.
        
        TODO: docstring"""

        coo = sparse.coo_matrix(self.X)
        data = OrderedDict([(0, coo.row+1), (1, coo.col+1), (2, coo.data)])
        df = pd.DataFrame(data, columns=data.keys())
        with open(file_path, 'w') as ofh:
            ofh.write('%%MatrixMarket matrix coordinate real general\n')
            ofh.write('%%%s\n' % '\t'.join(self.genes.astype(str)))
            ofh.write('%%%s\n' % '\t'.join(self.cells.astype(str)))
            ofh.write('%\n')
            ofh.write('%d %d %d\n' % (coo.shape[0], coo.shape[1], coo.nnz))
            df.to_csv(ofh, sep=' ', float_format='%.5f',
                      header=None, index=None)


    @classmethod
    def read_sparse(cls, file_path: str):
        """Read a sparse representation from a tab-delimited text file.
        
        TODO: docstring"""
        
        with open(file_path) as fh:
            next(fh)  # skip header line
            genes = next(fh)[1:-1].split('\t')
            cells = next(fh)[1:-1].split('\t')
            next(fh)
            m, n, nnz = [int(s) for s in next(fh)[:-1].split(' ')]
        
        t = pd.read_csv(file_path, sep=' ', skiprows=5, header=None,
                        dtype={0: np.uint32, 1: np.uint32})
        
        i = t[0].values - 1
        j = t[1].values - 1
        data = t[2].values

        assert data.size == nnz

        X = sparse.coo_matrix((data, (i,j)), shape=[m, n]).todense()

        return cls(X=X, genes=genes, cells=cells)

    @classmethod
    def read_10xgenomics(cls, tarball_fpath: str, prefix: str,
                         use_ensembl_ids: bool = False):
        """Read a 10X genomics compressed tarball containing expression data.
        
        Note: common prefix patterns:
        - "filtered_gene_bc_matrices/[annotations]/"
        - "filtered_matrices_mex/[annotations]/"

        TODO: docstring"""

        _LOGGER.info('Reading file: %s', tarball_fpath)

        with tarfile.open(tarball_fpath, mode='r:gz') as tf:
            ti = tf.getmember('%smatrix.mtx' % prefix)
            with tf.extractfile(ti) as fh:
                mtx = scipy.io.mmread(fh)

            ti = tf.getmember('%sgenes.tsv' % prefix)
            with tf.extractfile(ti) as fh:
                wrapper = io.TextIOWrapper(fh, encoding='ascii')
                i = 1
                if use_ensembl_ids:
                    i = 0
                gene_names = \
                        [row[i] for row in csv.reader(wrapper, delimiter='\t')]

            ti = tf.getmember('%sbarcodes.tsv' % prefix)
            with tf.extractfile(ti) as fh:
                wrapper = io.TextIOWrapper(fh, encoding='ascii')
                barcodes = \
                        [row[0] for row in csv.reader(wrapper, delimiter='\t')]


            assert mtx.shape[0] == len(gene_names)
            assert mtx.shape[1] == len(barcodes)
        
        _LOGGER.info('Matrix dimensions: %s', str(mtx.shape))
        X = mtx.todense()
        matrix = cls(X=X, genes=gene_names, cells=barcodes)
        
        return matrix
