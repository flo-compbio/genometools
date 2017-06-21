# Copyright (c) 2015-2017 Florian Wagner
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

"""Functions for configuring command-line arguments of GenomeTools scripts.
"""

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import *
import six

import sys
import argparse
# import textwrap

import genometools

file_mv = '<file>'
int_mv = '<int>'
float_mv = '<float>'
name_mv = '<name>'
str_mv = '<str>'

SYS_ENCODING = sys.getfilesystemencoding()


def str_type(a):
    """A string type for the `argparse.add_argument` ``type`` argument.

    This function fixes some issue where Python 2.7's argparse does not accept
    unicode values passed as arguments that are of the `future.types.newstr`
    type.
    """
    if six.PY2:
        return str(a, SYS_ENCODING)
    else:
        return str(a)


def get_argument_parser(prog=None, desc=None, formatter_class=None):
    """Create an argument parser.

    Parameters
    ----------
    prog: str
        The program name.
    desc: str
        The program description.
    formatter_class: argparse formatter class, optional
        The argparse formatter class to use.

    Returns
    -------
    `argparse.ArgumentParser`
        The arguemnt parser created.
    """
    if formatter_class is None:
        formatter_class = argparse.RawTextHelpFormatter

    parser = argparse.ArgumentParser(
        prog=prog, description=desc,
        formatter_class=formatter_class, add_help=False
    )

    g = parser.add_argument_group('Help')
    g.add_argument('-h', '--help', action='help',
                   help='Show this help message and exit.')

    v = genometools.__version__
    g.add_argument('--version', action='version', version='GenomeTools ' + v,
                   help='Output the GenomeTools version and exit.')

    return parser


def add_reporting_args(parser):
    """Add reporting arguments to an argument parser.

    Parameters
    ----------
    parser: `argparse.ArgumentParser`

    Returns
    -------
    `argparse.ArgumentGroup`
        The argument group created.
    """
    g = parser.add_argument_group('Reporting options')

    g.add_argument(
        '-l', '--log-file', default=None,
        type=str, metavar=file_mv,
        help='Path of log file (if specified, report to stdout AND file.'
    )

    g.add_argument('-q', '--quiet', action='store_true',
                   help='Only output errors and warnings.')

    g.add_argument(
        '-v', '--verbose', action='store_true',
        help='Enable verbose output. Ignored if --quiet is specified.'
    )

    return g
