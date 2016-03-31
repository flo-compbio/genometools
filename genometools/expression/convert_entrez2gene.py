#!/usr/bin/env python2.7

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

"""Script to Convert Entrez IDs to gene symbols in a gene expression matrix.

"""

import sys
import os
import textwrap
from collections import Counter, OrderedDict

import numpy as np

from genometools import misc
from genometools import cli
from genometools.expression import ExpGenome, ExpMatrix

def get_argument_parser():
    """Function to obtain the argument parser.

    Parameters
    ----------
    None

    Returns
    -------
    `argparse.ArgumentParser`
        A fully configured `argparse.ArgumentParser` object.

    Notes
    -----
    This function can also be used by the `sphinx-argparse` extension for
    sphinx to generate documentation for this script.
    """
    desc = 'Convert Entrez IDs to gene symbols.'
    parser = cli.get_argument_parser(desc = desc)

    file_mv = cli.file_mv

    g = parser.add_argument_group('Input and output files')

    g.add_argument('-e', '--expression-file', required = True,
            type = cli.str_type, metavar = file_mv,
            help = 'The expression file.')

    g.add_argument('-g', '--gene-file', required = True,
            type = cli.str_type, metavar = file_mv,
            help = textwrap.dedent('''\
                The gene file (e.g., generated by the
                ensembl_extract_protein_coding_genes.py script).'''))

            #help = 'The gene file (e.g,. generated by the ' +
            #       '`ensembl_extract_protein_coding_genes.py script).')

    g.add_argument('-c', '--entrez2gene-file', required = True,
            type = cli.str_type, metavar = file_mv,
            help = textwrap.dedent('''\
                The entrez2gene file (.e.g., generated by the
                ncbi_extract_entrez2gene.py script).'''))

    g.add_argument('-o', '--output-file', type = cli.str_type, required = True,
            metavar = file_mv, help = 'The output file.')

    g = parser.add_argument_group('Conversion options')

    g.add_argument('-s', '--strip-affy-suffix', action = 'store_true',
            help = textwrap.dedent('''\
                Strip the suffix "_at" from all Entrez IDs.
                (For use in affymetrix microarray pipeline.)'''))

    cli.add_reporting_args(parser)

    return parser

def main(args = None):

    if args is None:
        parser = get_argument_parser()
        args = parser.parse_args()

    expression_file = args.expression_file
    entrez2gene_file = args.entrez2gene_file
    gene_file = args.gene_file
    output_file = args.output_file

    strip_affy_suffix = args.strip_affy_suffix
    
    log_file = args.log_file
    quiet = args.quiet
    verbose = args.verbose

    # configure root logger
    logger = misc.get_logger(log_file = log_file, quiet = quiet,
            verbose = verbose)

    # read data
    genome = ExpGenome.read_tsv(gene_file)
    E = ExpMatrix.read_tsv(expression_file)
    e2g = dict(misc.read_all(entrez2gene_file))

    entrez = E.genes

    if strip_affy_suffix:
        # remove "_at" suffix from Entrez IDs
        entrez = [e[:-3] for e in entrez]
    logger.debug(str(entrez[:3]))

    # check that Entrez IDs are unique
    assert len(entrez) == len(set(entrez))

    # convert Entrez IDs to gene names
    f = 0
    genes = []
    X = []
    g = None
    for i,e in enumerate(entrez):
        #print e
        try:
            g = e2g[e]
        except KeyError:
            f += 1
        else:
            #assert g not in genes # check if there are multiple entrez IDs pointing to the same gene
            genes.append(g)
            X.append(E.X[i,:])
    assert len(genes) == len(set(genes))
    if f > 0:
        logger.warning('Failed to convert %d / %d entrez IDs ' +
                'to gene symbols (%.1f%%).',
                f, E.p, 100 * (f / float(E.p)))

    # filter for known protein-coding genes
    X = np.float64(X)
    p = X.shape[0]
    logger.debug(str(X.shape))
    sel = np.zeros(p, dtype=np.bool_)
    for i in range(p):
        if genes[i] in genome:
            sel[i] = True
    sel = np.nonzero(sel)[0]
    genes = [genes[i] for i in sel]
    X = X[sel,:]
    f = p - sel.size
    if f > 0:
        logger.warning('Failed to find %d / %d gene symbols in list of ' +
                'protein-coding genes (%.1f%%)',
                f, p, 100 * (f / float(p)))

    # generate new matrix (this automatically sorts the genes alphabetically)
    E_conv = ExpMatrix(genes, E.samples, X)

    # write output file
    E_conv.write_tsv(output_file)
 
    return 0

if __name__ == '__main__':
    return_code = main()
    sys.exit(return_code)
