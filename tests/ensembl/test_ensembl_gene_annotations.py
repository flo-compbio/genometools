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

"""Tests for functions in `cluster` module."""

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)

import os

import pytest

from genometools import ensembl
from genometools import misc

logger = misc.get_logger()

@pytest.mark.online
def test_latest_release():
    release = ensembl.get_latest_release()
    assert isinstance(release, int)
    logger.info('Current release: %d', release)

@pytest.mark.online
@pytest.mark.linux
@pytest.mark.darwin
@pytest.mark.cygwin
def test_download(my_download_dir):
    species = [
        'Homo_sapiens',
        'Mus_musculus',
    ]
    #dl_files = ensembl.download_gene_annotations(species, my_download_dir)
    species_data = ensembl.get_annotation_urls_and_checksums(species)
    assert len(species_data) == len(species)
    for spec, (url, chksum) in species_data.items():
        assert isinstance(url, str)
        assert isinstance(chksum, int)
    #for fn in dl_files:
    #    assert os.path.isfile(fn)