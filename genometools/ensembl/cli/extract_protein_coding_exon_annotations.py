#!/usr/bin/env python

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

"""Script for extracting all exons annotations of protein-coding genes.
"""

import sys
import os
import re
import argparse
import logging

import unicodecsv as csv

from genometools import misc
from genometools import ensembl
from genometools.ensembl.util import get_gtf_argument_parser
from genometools.gtf import parse_attributes

logger = logging.getLogger(__name__)


def get_argument_parser():
    """Function to obtain the argument parser.

    Returns
    -------
    A fully configured `argparse.ArgumentParser` object.

    Notes
    -----
    This function is used by the `sphinx-argparse` extension for sphinx.
    """
    desc = 'Extract all exon annotations of protein-coding genes.'
    parser = get_gtf_argument_parser(desc, default_field_name = 'exon')
    return parser

def main(args=None):
    """Extract all exon annotations of protein-coding genes."""

    if args is None:
        parser = get_argument_parser()
        args = parser.parse_args()

    input_file = args.annotation_file
    output_file = args.output_file
    species = args.species
    chrom_pat = args.chromosome_pattern
    field_name = args.field_name

    log_file = args.log_file
    quiet = args.quiet
    verbose = args.verbose

    # configure root logger
    log_stream = sys.stdout
    if output_file == '-':
        # if we print output to stdout, redirect log messages to stderr
        log_stream = sys.stderr

    logger = misc.get_logger(log_stream = log_stream, log_file = log_file,
            quiet = quiet, verbose = verbose)

    if chrom_pat is None:
        chrom_pat = re.compile(ensembl.SPECIES_CHROMPAT[species])
    else:
        chrom_pat = re.compile(chrom_pat)

    logger.info('Regular expression used for filtering chromosome names: "%s"',
            chrom_pat.pattern)

    chromosomes = set()
    excluded_chromosomes = set()
    i = 0
    exons = 0
    logger.info('Parsing data...')
    if input_file == '-':
        input_file = None
    with misc.smart_open_read(input_file, mode = 'rb', try_gzip = True) as fh, \
            misc.smart_open_write(output_file) as ofh:
        #if i >= 500000: break
        reader = csv.reader(fh, dialect = 'excel-tab')
        writer = csv.writer(ofh, dialect = 'excel-tab', lineterminator = os.linesep,
                quoting = csv.QUOTE_NONE , quotechar = '|')
        for l in reader:
            i += 1
            #if i % int(1e5) == 0:
            #   print '\r%d...' %(i), ; sys.stdout.flush() # report progress
            if len(l) > 1 and l[2] == field_name:
                attr = parse_attributes(l[8])
                type_ = attr['gene_biotype']
                if type_ in ['protein_coding','polymorphic_pseudogene']:

                    # test whether chromosome is valid
                    chrom = l[0]
                    m = chrom_pat.match(chrom)
                    if m is None:
                        excluded_chromosomes.add(chrom)
                        continue

                    chromosomes.add(chrom)
                    writer.writerow(l)
                    exons += 1

    logger.info('Done! (Parsed %d lines.)', i)

    logger.info('')
    logger.info('Gene chromosomes (%d):', len(chromosomes))
    logger.info('\t' + ', '.join(sorted(chromosomes)))
    logger.info('')
    logger.info('Excluded chromosomes (%d):', len(excluded_chromosomes))
    logger.info('\t' + ', '.join(sorted(excluded_chromosomes)))
    logger.info('')
    logger.info('Total no. of exons: %d' %(exons))

    return 0

if __name__ == '__main__':
    return_code = main()
    sys.exit(return_code)
