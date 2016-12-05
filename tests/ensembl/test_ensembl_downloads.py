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

"""Tests for functions that access the Enesembl FTP server."""

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import str as text
from builtins import int as newint

import os

import pytest

from genometools import ensembl

from genometools import misc

logger = misc.get_logger()

@pytest.mark.online
def test_latest_release():
    release = ensembl.get_latest_release()
    assert isinstance(release, newint)
    logger.info('Current release: %d', release)


@pytest.mark.online
def test_cdna_url():
    species = 'homo_sapiens' 
    url = ensembl.get_cdna_url(species, release=None, ftp=None)
    assert isinstance(url, text) and url != ''