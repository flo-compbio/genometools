Overview
========

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
