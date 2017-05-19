#!/usr/bin/env python

# Copyright (c) 2015, 2016 Florian Wagner
#
# This file is part of GenomeTools.
#
# GenomeTools is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License, Version 3,
# as published by the Free Software Foundation.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

"""Script for extracting lists of protein-coding genes.

This script parses an Ensembl GTF file containing _genome annotations, extracts
information about all protein-coding genes contained in it, and writes the
results to a tab-delimited text file. Each row in the output file corresponds
to one protein-coding gene.

Examples
--------

Extract the protein coding genes from the human Ensembl v82 gene annotations,
downloaded from the
`Ensembl FTP server <ftp://ftp.ensembl.org/pub/release-82/gtf/homo_sapiens/>`_:

.. code-block:: bash

    $ extract_protein_coding_genes.py \\
        -a Homo_sapiens.GRCh38.82.gtf.gz \\
        -o protein_coding_genes_human.tsv

"""

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
_oldstr = str
from builtins import *

import sys
import os
import re
import argparse
import logging
from collections import Counter

from ... import misc
from ... import ensembl
from ...gtf import parse_attributes
from .util import get_gtf_argument_parser


def get_argument_parser():
    """Function to obtain the argument parser.

    Returns
    -------
    A fully configured `argparse.ArgumentParser` object.

    Notes
    -----
    This function is used by the `sphinx-argparse` extension for sphinx.
    """
    desc = 'Extract all protein-coding genes from an Ensembl GTF file.'
    parser = get_gtf_argument_parser(desc)
    return parser


def main(args=None):
    """Extract protein-coding genes and store in tab-delimited text file.

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
    if not vinfo >= (2, 7):
        raise SystemError('Python interpreter version >= 2.7 required, '
                          'found %d.%d instead.' %(vinfo.major, vinfo.minor))

    if args is None:
        # parse command-line arguments
        parser = get_argument_parser()
        args = parser.parse_args()

    input_file = args.annotation_file
    output_file = args.output_file
    # species = args.species
    chrom_pat = args.chromosome_pattern
    log_file = args.log_file
    quiet = args.quiet
    verbose = args.verbose

    # configure root logger
    log_stream = sys.stdout
    if output_file == '-':
        # if we print output to stdout, redirect log messages to stderr
        log_stream = sys.stderr

    logger = misc.get_logger(log_stream=log_stream, log_file=log_file,
                             quiet=quiet, verbose=verbose)

    #if chrom_pat is None:
    #    chrom_pat = ensembl.SPECIES_CHROMPAT[species]

    if chrom_pat is not None:
        logger.info('Regular expression used for filtering chromosome names: '
                    '"%s"', chrom_pat)

    if input_file == '-':
        input_file = sys.stdin

    if output_file == '-':
        output_file = sys.stdout

    genes = ensembl.get_protein_coding_genes(
        input_file,
        chromosome_pattern=chrom_pat)
    genes.to_csv(output_file, sep='\t', index=False)

    return 0


if __name__ == '__main__':
    return_code = main()
    sys.exit(return_code)
