#!/usr/bin/env python2.7

# Copyright (c) 2015 Florian Wagner
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

"""Script for summarizing StringTie expression at the gene level.

"""

import sys
import os
import csv
#import gzip
import argparse
import logging

import numpy as np

from genometools import misc
from genometools.gtf import parse_attributes

def get_argument_parser():
    """Function to obtain the argument parser.

    Returns
    -------
    A fully configured `argparse.ArgumentParser` object.

    Notes
    -----
    This function is used by the `sphinx-argparse` extension for sphinx.

    """
    parser = argparse.ArgumentParser(description=
        'Extracts gene-level expression data from StringTie output.')

    parser.add_argument('-s','--stringtie-file', required=True,
        help="""Path of the StringTie output file .""")

    parser.add_argument('-g','--gene-file', required=True,
        help="""File containing a list of protein-coding genes.""")

    parser.add_argument('-o','--output-file', required=True,
        help="""Path of output file.""")

    parser.add_argument('-l','--log-file', default=None,
        help='Path of log file. If not specified, print to stdout.')

    parser.add_argument('-q','--quiet', action='store_true',
        help='Suppress all output except warnings and errors.')

    parser.add_argument('-v','--verbose', action='store_true',
        help='Enable verbose output. Ignored if ``--quiet`` is specified.')

    return parser

def main(args=None):
    """Extracts gene-level expression data from StringTie output.

    Parameters
    ----------
    args: argparse.Namespace object, optional
        The argument values. If not specified, the values will be obtained by
        parsing the command line arguments using the `argparse` module.

    Returns
    -------
    int
        Exit code (0 if no error occurred).
 
    """

    if args is None:
        # parse command-line arguments
        parser = get_argument_parser()
        args = parser.parse_args()

    stringtie_file = args.stringtie_file
    gene_file = args.gene_file
    output_file = args.output_file

    quiet = args.quiet
    verbose = args.verbose

    log_file = args.log_file
    log_stream = sys.stdout
    log_level = logging.INFO
    if quiet:
        log_level = logging.WARNING
    elif verbose:
        log_level = logging.DEBUG
    logger = misc.configure_logger(__name__, log_stream = log_stream,
            log_file = log_file, log_level = log_level)

    # read list of gene symbols
    genes = misc.read_single(gene_file)

    # read StringTie output file and summarize FPKM and TPM per gene
    n = len(genes)
    fpkm = np.zeros(n,dtype=np.float64)
    tpm = np.zeros(n,dtype=np.float64)
    with open(stringtie_file) as fh:
        reader = csv.reader(fh,dialect='excel-tab')
        for l in reader:
            if l[0][0] == '#':
                continue
            assert len(l) == 9
            if l[2] != 'transcript':
                continue
            attr = parse_attributes(l[8])
            g = attr['ref_gene_name']
            try:
                idx = misc.bisect_index(genes,g)
            except ValueError:
                logger.warning('Unknown gene: "%s".', g)
                continue
            f = float(attr['FPKM'])
            t = float(attr['TPM'])
            fpkm[idx] += f
            tpm[idx] += t

    # write output file
    E = np.c_[fpkm,tpm]
    with open(output_file,'w') as ofh:
        writer = csv.writer(ofh, dialect='excel-tab', lineterminator='\n',
                quoting=csv.QUOTE_NONE)
        for i,g in enumerate(genes):
            writer.writerow([g] + ['%.5f' %(e) for e in E[i,:]])

    return 0

if __name__ == '__main__':
    return_code = main()
    sys.exit(return_code)
