# Copyright (c) 2015, 2016 Florian Wagner
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

"""Utility functions for the Ensembl command-line interface."""

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
_oldstr = str
from builtins import *

import textwrap

from ... import ensembl
from ... import cli


def get_gtf_argument_parser(desc, default_field_name='gene'):
    """Return an argument parser with basic options for reading GTF files.

    Parameters
    ----------
    desc: str
        Description of the ArgumentParser
    default_field_name: str, optional
        Name of field in GTF file to look for.

    Returns
    -------
    `argparse.ArgumentParser` object
        The argument parser.
    """
    parser = cli.get_argument_parser(desc=desc)

    parser.add_argument(
        '-a', '--annotation-file', default='-', type=str,
        help=textwrap.dedent("""\
            Path of Ensembl gene annotation file (in GTF format). The file
            may be gzip'ed. If set to ``-``, read from ``stdin``.""")
    )

    parser.add_argument(
        '-o', '--output-file', required=True, type=str,
        help=textwrap.dedent("""\
            Path of output file. If set to ``-``, print to ``stdout``,
            and redirect logging messages to ``stderr``.""")
    )

    #parser.add_argument(
    #    '-s', '--species', type=str,
    #    choices=sorted(ensembl.SPECIES_CHROMPAT.keys()), default='human',
    #    help=textwrap.dedent("""\
    #        Species for which to extract genes. (This parameter is ignored
    #        if ``--chromosome-pattern`` is specified.)""")
    #)

    parser.add_argument(
        '-c', '--chromosome-pattern', type=str, required=False,
        default=None, help=textwrap.dedent("""\
            Regular expression that chromosome names have to match. [None]
            """)
    )

    #parser.add_argument(
    #    '-f', '--field-name', type=str, default=default_field_name,
    #    help=textwrap.dedent("""\
    #        Rows in the GTF file that do not contain this value
    #        in the third column are ignored.""")
    #)

    cli.add_reporting_args(parser)

    return parser
