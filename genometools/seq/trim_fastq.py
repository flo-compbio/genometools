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

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import *

import sys
# import argparse
# import itertools as it

from genometools import cli


def get_argument_parser():

    desc = 'Trim FASTQ file (read from stdin, write to stdout).'
    parser = cli.get_argument_parser(desc=desc)

    parser.add_argument(
        '-l', '--left', type=int, default=0,
        help="""Number of base pairs to trim from the left."""
    )

    parser.add_argument(
        '-r', '--right', type=int, default=0,
        help="""Number of base pairs to trim from the right."""
    )

    return parser


def main(args=None):

    if args is None:
        parser = get_argument_parser()
        args = parser.parse_args()

    left = args.left
    right = args.right

    # counter = it.cycle(range(4))
    from_right = right + 1
    while True:
        try:
            sys.stdout.write(next(sys.stdin))
        except StopIteration:
            break
        sys.stdout.write(next(sys.stdin)[left:-from_right] + '\n')
        sys.stdout.write(next(sys.stdin))
        sys.stdout.write(next(sys.stdin)[left:-from_right] + '\n')

    return 0

if __name__ == '__main__':
    return_code = main()
    sys.exit(return_code)
