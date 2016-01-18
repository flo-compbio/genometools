.. GenomeTools documentation master file, created by
   sphinx-quickstart on Fri Sep 11 14:49:21 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

GenomeTools |version|
=====================

GenomeTools is a collection of Python 2.7 scripts and functions that perform
general tasks encountered in the analysis of genomic data. GenomeTools is free
and open-source software (see `License <license>`).

Main Features
-------------

- Script to extract a list of all protein-coding genes from Ensembl GTF files: :ref:`extract_genes`
- Script to extract a mapping from Entrez IDs to gene symbols: :ref:`extract_e2g`
- Basic classes for working with genomic data: `genometools.basic <basic>`
- Classes for working with expression data: `genometools.expression <expression>`
- Miscellaneous convenience functions: :doc:`genometools.misc <misc>`

Demo Notebooks
--------------

- `Scripts.ipynb`__ (:download:`download <notebooks/Scripts.ipynb>`): GenomeTools command-line scripts

__ scripts_notebook_

.. _scripts_notebook: http://nbviewer.ipython.org/url/genometools.readthedocs.org/en/latest/_downloads/Scripts.ipynb

.. This only links to the "latest" (master) version...no way to automatically switch to the "develop" version

Installation
------------

GenomeTools can be installed from `PyPI <https://pypi.python.org/pypi>`_ using `pip <https://pip.pypa.io/en/stable/>`_:

.. code-block:: bash

    $ pip install genometools

.. A link to :mod:`genometools.misc.misc` and another link to :mod:`misc` and one to `misc`.
.. A link to :mod:`genometools.ensembl.extract_protein_coding_genes` and another link to :mod:`misc` and one to `extract_protein_coding_genes`.

.. Indices and tables
.. ==================
.. 
.. * :ref:`genindex`
.. * :ref:`search`
.. * :ref:`modindex`

.. toctree::
    :maxdepth: 2
    :hidden:

    self
    Scripts <scripts>
    Basic classes <basic>
    Expression data <expression>
    Miscellaneous <misc>
    license
