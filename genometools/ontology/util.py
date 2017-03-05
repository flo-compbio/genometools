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

"""Functions for downloading GO annotations"""

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
_oldstr = str
from builtins import *

# import os
# import sys
# import ftplib
# import re
# import tempfile
# from collections import Iterable
from contextlib import closing

import pandas as pd
import requests

from .. import misc

logger = misc.get_logger()

#ftp_server = 'ftp.ensembl.org'
#user = 'anonymous'


def get_latest_release():
    """Gets the name (date) of the latest release of the Gene Ontology."""
    list_url = 'http://viewvc.geneontology.org/viewvc/GO-SVN/ontology-releases/'
    with closing(requests.get(list_url)) as r:
        text = r.text
    all_versions = re.findall('<a name="(\d{4}-\d\d-\d\d)" href="', text)
    latest = list(sorted(all_versions))[-1]
    return latest


def download_release(download_file, release=None):
    """Downloads the "go-basic.obo" file for the specified release."""
    if release is None:
        release = get_latest_release()
    url = 'http://viewvc.geneontology.org/viewvc/GO-SVN/ontology-releases/%s/go-basic.obo' % release
    #download_file = 'go-basic_%s.obo' % release
    misc.http_download(url, download_file)


def get_current_ontology_date():
    """Get the release date of the current Gene Ontolgo release."""
    with closing(requests.get(
            'http://geneontology.org/ontology/go-basic.obo',
            stream=True)) as r:
        for i, l in enumerate(r.iter_lines(decode_unicode=True)):
            if i == 1:
                assert l.split(':')[0] == 'data-version'
                date = l.split('/')[-1]
                break

    return date
