#!/usr/bin/env python

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

"""Script to extract a mapping of Entrez IDs to gene symbols.

This script parses the gene2accession.gz file from the `NCBI FTP server`__
(or a filtered version thereof), extracts a mapping of Entrez IDs to gene
symbols, and writes this mapping to a tab-delimited text file. Each row in the
output file contains one Entrez ID and its associated gene symbol.

__ ncbi_ftp_

Examples
--------

- Step 1) Filtering of the gene2accession.gz file to only contain entries for
  human genes:

.. code-block:: bash

    $ gunzip -c gene2accession.gz | grep -P "^9606\t" | \
            gzip > gene2accession_human.gz

- Step 2) Use GenomeTools to extract a mapping between Entrez IDs and gene
  symbols:

.. code-block:: bash

    $ extract_entrez2gene.py \\
        -f gene2accession_human.gz \\
        -o entrez2gene_human.tsv

.. _ncbi_ftp: ftp://ftp.ncbi.nlm.nih.gov/gene/DATA

"""

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import *

import sys
import os
# import logging
# import argparse
# import gzip
import textwrap

import unicodecsv as csv

from ... import misc
from ... import cli


def get_argument_parser():
    """Creates the argument parser for the extract_entrez2gene.py script.

    Returns
    -------
    A fully configured `argparse.ArgumentParser` object.

    Notes
    -----
    This function is used by the `sphinx-argparse` extension for sphinx.

    """

    desc = 'Generate a mapping of Entrez IDs to gene symbols.'

    parser = cli.get_argument_parser(desc=desc)

    parser.add_argument(
        '-f', '--gene2acc-file', type=str, required=True,
        help=textwrap.dedent("""\
            Path of gene2accession.gz file (from
            ftp://ftp.ncbi.nlm.nih.gov/gene/DATA), or a filtered version
            thereof.""")
    )

    parser.add_argument(
        '-o', '--output-file', type=str, required=True,
        help=textwrap.dedent("""\
            Path of output file. If set to ``-``, print to ``stdout``,
            and redirect logging messages to ``stderr``.""")
    )

    parser.add_argument(
        '-l', '--log-file', type=str, default=None,
        help='Path of log file. If not specified, print to stdout.'
    )

    parser.add_argument(
        '-q', '--quiet', action='store_true',
        help='Suppress all output except warnings and errors.'
    )

    parser.add_argument(
        '-v', '--verbose', action='store_true',
        help='Enable verbose output. Ignored if ``--quiet`` is specified.'
    )

    return parser


def read_gene2acc(file_path, logger):
    """Extracts Entrez ID -> gene symbol mapping from gene2accession.gz file.

    Parameters
    ----------
    file_path: str
        The path of the gene2accession.gz file (or a filtered version thereof).
        The file may be gzip'ed.

    Returns
    -------
    dict
        A mapping of Entrez IDs to gene symbols.
    """
    gene2acc = {}
    with misc.smart_open_read(file_path, mode='rb', try_gzip=True) as fh:
        reader = csv.reader(fh, dialect='excel-tab')
        next(reader)  # skip header
        for i, l in enumerate(reader):
            id_ = int(l[1])
            symbol = l[15]

            try:
                gene2acc[id_].append(symbol)
            except KeyError:
                gene2acc[id_] = [symbol]

            # print (l[0],l[15])

    # make sure all EntrezIDs map to a unique gene symbol
    n = len(gene2acc)
    for k, v in gene2acc.items():
        symbols = sorted(set(v))
        assert len(symbols) == 1
        gene2acc[k] = symbols[0]

    all_symbols = sorted(set(gene2acc.values()))
    m = len(all_symbols)

    logger.info('Found %d Entrez Gene IDs associated with %d gene symbols.',
                n, m)
    return gene2acc


def write_entrez2gene(file_path, entrez2gene, logger):
    """Writes Entrez ID -> gene symbol mapping to a tab-delimited text file.

    Parameters
    ----------
    file_path: str
        The path of the output file.
    entrez2gene: dict
        The mapping of Entrez IDs to gene symbols.

    Returns
    -------
    None

    """
    with misc.smart_open_write(file_path, mode='wb') as ofh:
        writer = csv.writer(ofh, dialect='excel-tab',
                            lineterminator=os.linesep)
        for k in sorted(entrez2gene.keys(), key=lambda x: int(x)):
            writer.writerow([k, entrez2gene[k]])
    logger.info('Output written to file "%s".', file_path)


def main(args=None):
    """Extracts Entrez ID -> gene symbol mapping and writes it to a text file.

    Parameters
    ----------
    args: argparse.Namespace object, optional
        The argument values. If not specified, the values will be obtained by
        parsing the command line arguments using the `argparse` module.

    Returns
    -------
    int
        Exit code (0 if no error occurred).

    Raises
    ------
    SystemError
        If the version of the Python interpreter is not >= 2.7.
    """
    vinfo = sys.version_info
    if not (vinfo >= (2, 7)):
        raise SystemError('Python interpreter version >= 2.7 required, '
                          'found %d.%d instead.' % (vinfo.major, vinfo.minor))

    if args is None:
        parser = get_argument_parser()
        args = parser.parse_args()

    gene2acc_file = args.gene2acc_file
    output_file = args.output_file
    log_file = args.log_file
    quiet = args.quiet
    verbose = args.verbose

    # configure logger
    log_stream = sys.stdout
    if output_file == '-':
        log_stream = sys.stderr

    logger = misc.get_logger(log_stream=log_stream, log_file=log_file,
                             quiet=quiet, verbose=verbose)

    entrez2gene = read_gene2acc(gene2acc_file, logger)
    write_entrez2gene(output_file, entrez2gene, logger)

    return 0

if __name__ == '__main__':
    return_code = main()
    sys.exit(return_code)
