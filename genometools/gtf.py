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

"""Functions for GTF files."""

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import *

import re

def parse_attributes(s):
    """ Parses the ``attribute`` string of a GFF/GTF annotation.

    Parameters
    ----------
    s : str
        The attribute string.

    Returns
    -------
    dict
        A dictionary containing attribute name/value pairs.

    Notes
    -----
    The ``attribute`` string is the 9th field of each annotation (row),
    as described in the
    `GTF format specification <http://mblab.wustl.edu/GTF22.html>`_.
    """
    # use regular expression with negative lookbehind to make sure we don't
    # split on escaped semicolons ("\;")
    attr_sep = re.compile(r'(?<!\\)\s*;\s*')
    attr = {}
    atts = attr_sep.split(s)
    for a in atts:
        #print(a)
        kv = a.split(' ', maxsplit=1)
        if len(kv) == 2:
            k, v = kv
            v = v.strip('"')
            attr[k] = v
    return attr 
