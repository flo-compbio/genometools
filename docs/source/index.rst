.. GenomeTools documentation master file, created by
   sphinx-quickstart on Fri Sep 11 14:49:21 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

GenomeTools |version|
=====================

GenomeTools is a collection of scripts and functions that perform general
tasks encountered in the analysis of genomic data. GenomeTools is free
and open-source software, `licensed <license>` under the GNU GPL v3.

Installation
------------

GenomeTools can be installed from `PyPI <https://pypi.python.org/pypi>`_ using `pip <https://pip.pypa.io/en/stable/>`_:

.. code-block:: bash

    $ pip install genometools

Main Features
-------------

- Script to extract a list of all protein-coding genes from Ensembl GTF files: :ref:`extract_genes`
- Script to extract a mapping from Entrez IDs to gene symbols: :ref:`extract_e2g`
- Convenience function for configuring a logger: `misc.get_logger`

.. A link to :mod:`genometools.misc` and another link to :mod:`misc` and one to `misc`.
.. A link to :mod:`genometools.extract_protein_coding_genes` and another link to :mod:`misc` and one to `extract_protein_coding_genes`.

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
    Modules <modules>
    Scripts <scripts>
    license
