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

"""Script for filtering FASTA file by chromosome names.

FASTA files are used to store genomic sequences and can be downloaded from
the Ensembl `download site`__. Ensembl formats them to contain 60 nucleotides
per line, and uses Unix-style line breaks (consisting of a single ``LF`` byte;
ASCII value 10 / 0x0A). The chromosome names are followed by a space and a
string that contains additional information about the chromosome.

__ ensembl_download_

.. _ensembl_download: http://www.ensembl.org/info/data/ftp/index.html

"""

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import *

import sys
# import os
import re
import textwrap
# import logging

from genometools import misc
from genometools import cli
from genometools import ensembl
from genometools.seq import FastaReader


def get_argument_parser():
    """Returns an argument parser object for the script."""

    desc = 'Filter FASTA file by chromosome names.'
    parser = cli.get_argument_parser(desc=desc)

    parser.add_argument(
        '-f', '--fasta-file', default='-', type=str, help=textwrap.dedent("""\
                Path of the FASTA file. The file may be gzip'ed.
                If set to ``-``, read from ``stdin``."""))

    parser.add_argument(
        '-s', '--species', type=str,
        choices=sorted(ensembl.SPECIES_CHROMPAT.keys()),
        default='human', help=textwrap.dedent("""\
            Species for which to extract genes. (This parameter is ignored
            if ``--chromosome-pattern`` is specified.)""")
    )

    parser.add_argument(
        '-c', '--chromosome-pattern', type=str, required=False,
        default=None, help=textwrap.dedent("""\
            Regular expression that chromosome names have to match.
            If not specified, determine pattern based on the setting of
            ``--species``.""")
    )

    parser.add_argument(
        '-o', '--output-file', type=str, required=True,
        help=textwrap.dedent("""\
            Path of output file. If set to ``-``, print to ``stdout``,
            and redirect logging messages to ``stderr``."""))

    parser = cli.add_reporting_args(parser)
    
    return parser


def main(args=None):
    """Script body."""

    if args is None:
        # parse command-line arguments 
        parser = get_argument_parser()
        args = parser.parse_args()

    fasta_file = args.fasta_file
    species = args.species
    chrom_pat = args.chromosome_pattern
    output_file = args.output_file
    
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

    # generate regular expression object from the chromosome pattern
    if chrom_pat is None:
        chrom_pat = ensembl.SPECIES_CHROMPAT[species]
    chrom_re = re.compile(chrom_pat)

    # filter the FASTA file
    # note: each chromosome sequence is temporarily read into memory,
    # so this script has a large memory footprint
    with \
        misc.smart_open_read(
            fasta_file, mode='r', encoding='ascii', try_gzip=True
        ) as fh, \
        misc.smart_open_write(
            output_file, mode='w', encoding='ascii'
        ) as ofh:

        # inside = False
        reader = FastaReader(fh)
        for seq in reader:
            chrom = seq.name.split(' ', 1)[0]
            if chrom_re.match(chrom) is None:
                logger.info('Ignoring chromosome "%s"...', chrom)
                continue
            seq.name = chrom
            seq.append_fasta(ofh)

    return 0

if __name__ == '__main__':
    return_code = main()
    sys.exit(return_code)
