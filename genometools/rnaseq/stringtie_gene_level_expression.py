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

"""Script for summarizing StringTie expression at the gene level.

"""

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import *

import sys
import os
# import argparse
# import logging
# from collections import Counter

import numpy as np

import unicodecsv as csv

from genometools import misc
from genometools import cli
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
    file_mv = cli.file_mv

    desc = 'Extracts gene-level expression data from StringTie output.'
    parser = cli.get_argument_parser(desc)

    parser.add_argument(
        '-s', '--stringtie-file', type=str, required=True, metavar=file_mv,
        help="""Path of the StringTie output file ."""
    )

    parser.add_argument(
        '-g', '--gene-file', type=str, required=True, metavar=file_mv,
        help="""File containing a list of protein-coding genes."""
    )

    parser.add_argument(
        '--no-novel-transcripts', action='store_true',
        help="""Ignore novel transcripts."""
    )

    # parser.add_argument(
    #     '--ambiguous-transcripts', default = 'ignore',
    #      help='Strategy for counting expression of ambiguous novel '
    #            'transcripts.'
    # )
    # possible strategies for ambiguous transcripts: 'ignore','highest','all'

    parser.add_argument(
        '-o', '--output-file', type=str, required=True, metavar=file_mv,
        help="""Path of output file."""
    )

    cli.add_reporting_args(parser)

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
    no_novel_transcripts = args.no_novel_transcripts
    output_file = args.output_file

    log_file = args.log_file
    quiet = args.quiet
    verbose = args.verbose

    logger = misc.get_logger(log_file=log_file, quiet=quiet, verbose=verbose)

    # read list of gene symbols
    logger.info('Reading gene data...')
    genes = misc.read_single(gene_file)

    # read StringTie output file and summarize FPKM and TPM per gene
    logger.info('Parsing StringTie output...')

    logger.info('Associating StringTie gene IDs with gene symbols...')
    stringtie_genes = {}
    with open(stringtie_file) as fh:
        reader = csv.reader(fh, dialect='excel-tab')
        for l in reader:
            if l[0][0] == '#':
                continue
            assert len(l) == 9
            if l[2] != 'transcript':
                continue
            attr = parse_attributes(l[8])
            try:
                ref_gene = attr['ref_gene_name']
            except KeyError:
                continue
            else:
                # entry has a "ref_gene_name" attribute
                try:
                    g = stringtie_genes[attr['gene_id']]
                except KeyError:
                    stringtie_genes[attr['gene_id']] = {ref_gene, }
                else:
                    g.add(ref_gene)
    logger.info('Associated %d gene IDs with gene symbols.',
                len(stringtie_genes))
    # C = Counter(len(v) for v in stringtie_genes.itervalues())
    gene_ids_ambiguous = [k for k, v in stringtie_genes.items()
                          if len(v) > 1]
    n = len(gene_ids_ambiguous)
    logger.info('%d / %d associated with multiple gene symbols (%.1f%%).',
                n, len(stringtie_genes), 100*(n/float(len(stringtie_genes))))

    # read StringTie output file and summarize FPKM and TPM per gene
    n = len(genes)
    fpkm = np.zeros(n, dtype=np.float64)
    tpm = np.zeros(n, dtype=np.float64)
    fpkm_novel_gene = 0
    fpkm_unknown_gene_name = 0
    fpkm_novel_trans = 0
    fpkm_ambig = 0
    with open(stringtie_file) as fh:
        reader = csv.reader(fh, dialect='excel-tab')
        for l in reader:
            if l[0][0] == '#':
                # skip header
                continue
            assert len(l) == 9

            if l[2] != 'transcript':
                # skip exon lines
                continue

            attr = parse_attributes(l[8])
            f = float(attr['FPKM'])

            try:
                g = attr['ref_gene_name']
            except KeyError:
                if no_novel_transcripts:
                    # ignore this transcript
                    fpkm_novel_trans += f
                    continue
                else:
                    # see if we can assign a gene name based on the gene ID
                    try:
                        assoc = stringtie_genes[attr['gene_id']]
                    except KeyError:
                        # gene_id not associated with any reference gene
                        fpkm_novel_gene += f
                        continue
                    else:
                        if len(assoc) > 1:
                            # gene ID associated with multiple ref. genes
                            # => ingored
                            fpkm_ambig += f
                            continue
                        else:
                            # gene ID associated with exactly one ref. gene
                            g = list(assoc)[0]
 
            try:
                idx = misc.bisect_index(genes, g)
            except ValueError:
                fpkm_unknown_gene_name += f
                logger.warning('Unknown gene name: "%s".', g)
                continue

            t = float(attr['TPM'])
            fpkm[idx] += f
            tpm[idx] += t

    # ignored_fpkm = None
    if no_novel_transcripts:
        ignored_fpkm = fpkm_novel_trans + fpkm_unknown_gene_name
    else:
        ignored_fpkm = fpkm_novel_gene + fpkm_ambig + fpkm_unknown_gene_name
    total_fpkm = np.sum(fpkm) + ignored_fpkm
    logger.info('Ignored %.1f / %.1f FPKM (%.1f%%)', ignored_fpkm,
                total_fpkm, 100*(ignored_fpkm/total_fpkm))

    if no_novel_transcripts and fpkm_novel_trans > 0:
        logger.info('Ignored %.1f FPKM from novel transcripts (%.1f%%).',
                    fpkm_novel_trans, 100*(fpkm_novel_trans/total_fpkm))

    else:
        if fpkm_novel_gene > 0:
            logger.info('Ignored %.1f FPKM from transcripts of novel genes '
                        '(%.1f%%).',
                        fpkm_novel_gene, 100*(fpkm_novel_gene/total_fpkm))

        if fpkm_ambig > 0:
            logger.info('Ignored %.1f FPKM from transcripts with ambiguous '
                        'gene membership (%.1f%%).',
                        fpkm_ambig, 100*(fpkm_ambig/total_fpkm))

    if fpkm_unknown_gene_name > 0:
        logger.info('Ignored %.1f FPKM from transcripts of genes with unknown '
                    'names (%.1f%%).',
                    fpkm_unknown_gene_name,
                    100*(fpkm_unknown_gene_name/total_fpkm))

    # write output file
    E = np.c_[fpkm, tpm]
    with open(output_file, 'w') as ofh:
        writer = csv.writer(ofh, dialect='excel-tab',
                            lineterminator=os.linesep,
                            quoting=csv.QUOTE_NONE)
        for i, g in enumerate(genes):
            writer.writerow([g] + ['%.5f' % e for e in E[i, :]])

    return 0

if __name__ == '__main__':
    return_code = main()
    sys.exit(return_code)
