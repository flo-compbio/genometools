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

"""Script to find all SRA run IDs (SRR*) for a given SRA experiment ID (SRX*).
"""

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import *

import sys
import os
# import argparse
# import logging
import ftputil

import unicodecsv as csv

from ... import misc
from ... import cli


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

    desc = 'Find all runs (SRR..) associated with an SRA experiment (SRX...).'

    parser = cli.get_argument_parser(desc=desc)

    parser.add_argument(
        '-e', '--experiment-file', type=str, required=True, metavar=file_mv,
        help='File with SRA experiment IDs (starting with "SRX").'
    )

    parser.add_argument(
        '-o', '--output-file', type=str, required=True, metavar=file_mv,
        help='The output file.'
    )

    cli.add_reporting_args(parser)

    return parser


def sizeof_fmt(num, suffix='B'):
    for unit in ['', 'Ki', 'Mi', 'Gi', 'Ti', 'Pi', 'Ei', 'Zi']:
        if abs(num) < 1024.0:
            return "%3.0f%s%s" % (num, unit, suffix)
        num /= 1024.0
    return "%.0f%s%s" % (num, 'Yi', suffix)


def main(args=None):
    """Download all .sra from NCBI SRA for a given experiment ID.

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

    experiment_file = args.experiment_file
    output_file = args.output_file

    # log_file = args.log_file
    # quiet = args.quiet
    # verbose = args.verbose

    # logger = misc.get_logger(log_file=log_file, quiet=quiet,
    #                          verbose=verbose)

    host = 'ftp-trace.ncbi.nlm.nih.gov'
    user = 'anonymous'
    password = 'anonymous'

    # output_dir = download_dir + experiment_id + '/'
    # make sure output directory exists
    # misc.make_sure_dir_exists(output_dir)
    # logger.info('Created output directory: "%s".', output_dir)

    experiments = misc.read_single(experiment_file)

    runs = []
    with ftputil.FTPHost(host, user, password) as ftp_host:
        for exp in experiments:
            exp_dir = '/sra/sra-instant/reads/ByExp/sra/SRX/%s/%s/' \
                    % (exp[:6], exp)
            ftp_host.chdir(exp_dir)
            run_folders = ftp_host.listdir(ftp_host.curdir)
            # logging.info('Found %d run folders.',len(run_folders))

            for folder in run_folders:
                files = ftp_host.listdir(folder)
                assert len(files) == 1
                runs.append((exp, folder))

    with open(output_file, 'wb') as ofh:
        writer = csv.writer(ofh, dialect='excel-tab',
                            lineterminator=os.linesep,
                            quoting=csv.QUOTE_NONE)
        for r in runs:
            writer.writerow(r)
        
    return 0

if __name__ == '__main__':
    return_code = main()
    sys.exit(return_code)
