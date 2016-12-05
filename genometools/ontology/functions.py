# Copyright (c) 2016 Florian Wagner
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

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
_oldstr = str
from builtins import *

"""Functions for downloading GO annotations"""

import os
import sys
import ftplib
import re
import tempfile
from collections import Iterable
from contextlib import closing

import pandas as pd
import requests

from .. import misc

logger = misc.get_logger()

#ftp_server = 'ftp.ensembl.org'
#user = 'anonymous'


def get_current_ontology_date():
        with closing(requests.get(
                'http://geneontology.org/ontology/go-basic.obo',
                stream=True)) as r:
            for i, l in enumerate(r.iter_lines(decode_unicode=True)):
                if i == 1:
                    assert l.split(':')[0] == 'data-version'
                    date = l.split('/')[-1]
                    break
        
        return date