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

"""Command-line interface to extract a mapping of Entrez IDs to gene symbols.

The script contained in the `main` function reads the gene2accession.gz file
from ftp://ftp.ncbi.nlm.nih.gov/gene/DATA (or a filtered version thereof) and
extracts a mapping of Entrez IDs to gene symbols.

Examples
--------

- Step 1) Filtering of the gene2accession.gz file to only contain entries for
  human genes:

.. code-block:: bash

    $ gunzip -c gene2accession.gz | grep -P "^9606\t" | gzip > gene2accession_human.gz

- Step 2) Use GenomeTools to extract a mapping between Entrez IDs and gene
  symbols:

.. code-block:: bash

    $ extract_entrez2gene.py \\
        -f gene2accession_human.gz \\
        -o entrez2gene_human.tsv

"""

import sys
import os
import logging
import argparse
import csv
import gzip

def get_argument_parser():
    """Creates the argument parser for the extract_entrez2gene.py script.

    Returns
    -------
    A fully configured `argparse.ArgumentParser` object.

    Notes
    -----
    This function is used by the `sphinx-argparse` extension for sphinx.

    """

    parser = argparse.ArgumentParser(
            description = 'Generate a mapping of Entrez IDs to gene symbols.')

    parser.add_argument('-f','--gene2acc-file',required=True,\
            help='Path of gene2accession.gz file \
            (from ftp://ftp.ncbi.nlm.nih.gov/gene/DATA),\
            or a filtered version thereof.')
    parser.add_argument('-o','--output-file',required=True,\
            help='Path of output file.')
    parser.add_argument('-l','--log-file',default=None,\
            help='Path of log file. If not specified, print to stdout.')
    parser.add_argument('-v','--verbose',action='store_true',\
            help='Enable verbose output.')

    return parser

def read_gene2acc(fn,logger):
    """Extracts Entrez ID -> gene symbol mapping from gene2accession.gz file.

    Parameters
    ----------
    fn: str
        The path to the gene2accession.gz file.
    logger: logging.Logger object
        The logger.

    Returns
    -------
    dict
        A mapping of Entrez IDs to gene symbols.

    """

    gene2acc = {}
    with open_plain_or_gzip(fn) as fh:
        reader = csv.reader(fh,dialect='excel-tab')
        reader.next() # skip header
        for i,l in enumerate(reader):
            id_ = int(l[1])
            symbol = l[15]

            try:
                gene2acc[id_].append(symbol)
            except:
                gene2acc[id_] = [symbol]
            #print l[0],l[15]

    # make sure all EntrezIDs map to a unique gene symbol
    n = len(gene2acc.keys())
    for k,v in gene2acc.iteritems():
        symbols = sorted(set(v))
        assert len(symbols) == 1
        gene2acc[k] = symbols[0]

    all_symbols = sorted(set(gene2acc.values()))
    m = len(all_symbols)

    logger.info('Found %d Entrez Gene IDs associated with %d gene symbols.', n,m)
    return gene2acc

def write_entrez2gene(ofn,entrez2gene,logger):
    """Writes Entrez ID -> gene symbol mapping to a tab-delimited text file.

    Parameters
    ----------
    ofn: str
        The path to the output file file.
    entrez2gene: dict
        The mapping of Entrez IDs to gene symbols.
    logger: logging.Logger object
        The logger.

    Returns
    -------
    None

    """

    with open(ofn,'w') as ofh:
        writer = csv.writer(ofh,dialect='excel-tab',lineterminator='\n')
        for k,v in entrez2gene.iteritems():
            writer.writerow([k,v])
    logger.info('Output written to file "%s".', ofn)

def main(args=None):
    """Extracts Entrez ID -> gene symbol mapping and writes it to a text file.

    This is the main function of the extract_entrez2gene.py script, which
    parses a gene2accession.gz file, extracts a mapping of Entrez IDs to gene
    symbols, and writes this mapping to a tab-delimited text file. Each row in
    the output file contains one Entrez ID ands its associated gene symbol.

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
        parser = get_argument_parser()
        args = parser.parse_args()

    gene2acc_file = args.gene2acc_file
    output_file = args.output_file
    log_file = args.log_file
    verbose = args.verbose

    # configure logger
    log_level = logging.INFO
    if verbose:
        log_level = logging.DEBUG

    log_format = '[%(asctime)s] %(levelname)s: %(message)s'
    log_datefmt = '%Y-%m-%d %H:%M:%S'
    # when filename is not None, "stream" parameter is ignored (see https://docs.python.org/2/library/logging.html#logging.basicConfig)
    logging.basicConfig(filename=log_file,stream=sys.stdout,level=log_level,format=log_format,datefmt=log_datefmt)
    logger = logging.getLogger()

    entrez2gene = read_gene2acc(gene2acc_file,logger)
    write_entrez2gene(output_file,entrez2gene,logger)

    return 0

if __name__ == '__main__':
    return_code = main()
    sys.exit(return_code)
