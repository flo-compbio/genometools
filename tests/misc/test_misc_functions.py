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

"""Tests for functions in `misc` module."""

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import str as text

import os

import pytest

from genometools import misc

logger = misc.get_logger()

@pytest.mark.linux
@pytest.mark.darwin
def test_checksum(my_checksum_file):
    assert misc.get_file_checksum(my_checksum_file) == 2761
    assert misc.test_file_checksum(my_checksum_file, 2761)


def test_ftp_download(my_readme_file):
    misc.ftp_download('ftp://ftp.ensembl.org/pub/current_README',                                 my_readme_file)
    assert os.stat(my_readme_file).st_size > 0