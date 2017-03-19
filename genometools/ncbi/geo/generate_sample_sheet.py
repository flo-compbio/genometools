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

"""Script for generating a sample sheet based on a GEO series matrix.

"""

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import *
import six

import sys
import os

import unicodecsv as csv

from ... import misc
# from genometools import ensembl
from ... import cli


def get_argument_parser():
    """Create the argument parser for the script.

    Parameters
    ----------

    Returns
    -------
    `argparse.ArgumentParser`
        The arguemnt parser.
    """
    desc = 'Generate a sample sheet based on a GEO series matrix.'
    parser = cli.get_argument_parser(desc=desc)

    g = parser.add_argument_group('Input and output files')

    g.add_argument(
        '-s', '--series-matrix-file', type=cli.str_type, required=True,
        metavar=cli.file_mv, help='The GEO series matrix file.'
    )

    g.add_argument(
        '-o', '--output-file', type=cli.str_type,
        required=True,
        metavar=cli.file_mv, help='The output file.'
    )

    g.add_argument(
        '-e', '--encoding', type=cli.str_type,
        metavar=cli.str_mv, default='UTF-8',
        help='The encoding of the series matrix file. [UTF-8]'
    )

    cli.add_reporting_args(parser)

    return parser


def read_series_matrix(path, encoding):
    """Read the series matrix."""
    assert isinstance(path, str)

    accessions = None
    titles = None
    celfile_urls = None
    with misc.smart_open_read(path, mode='rb', try_gzip=True) as fh:
        reader = csv.reader(fh, dialect='excel-tab', encoding=encoding)
        for l in reader:
            if not l:
                continue
            if l[0] == '!Sample_geo_accession':
                accessions = l[1:]
            elif l[0] == '!Sample_title':
                titles = l[1:]
            elif l[0] == '!Sample_supplementary_file' and celfile_urls is None:
                celfile_urls = l[1:]
            elif l[0] == '!series_matrix_table_begin':
                # we've read the end of the section containing metadata
                break
    return accessions, titles, celfile_urls


def write_sample_sheet(path, accessions, names, celfile_urls, sel=None):
    """Write the sample sheet."""
    with open(path, 'wb') as ofh:
        writer = csv.writer(ofh, dialect='excel-tab',
                            lineterminator=os.linesep,
                            quoting=csv.QUOTE_NONE)
        # write header
        writer.writerow(['Accession', 'Name', 'CEL file', 'CEL file URL'])
        n = len(names)
        if sel is None:
            sel = range(n)
        for i in sel:
            cf = celfile_urls[i].split('/')[-1]
            # row = [accessions[i], names[i], cf, celfile_urls[i]]
            writer.writerow([accessions[i], names[i], cf, celfile_urls[i]])


def main(args=None):
    """Script entry point."""

    if args is None:
        parser = get_argument_parser()
        args = parser.parse_args()

    #series_matrix_file = newstr(args.series_matrix_file, 'utf-8')
    #output_file = newstr(args.output_file, 'utf-8')
    #encoding = newstr(args.encoding, 'utf-8')
    series_matrix_file = args.series_matrix_file
    output_file = args.output_file
    encoding = args.encoding

    # log_file = args.log_file
    # quiet = args.quiet
    # verbose = args.verbose

    # logger = misc.get_logger(log_file = log_file, quiet = quiet,
    #        verbose = verbose)

    accessions, titles, celfile_urls = read_series_matrix(
        series_matrix_file, encoding=encoding)
    write_sample_sheet(output_file, accessions, titles, celfile_urls)

    return 0

if __name__ == '__main__':
    return_code = main()
    sys.exit(return_code)
