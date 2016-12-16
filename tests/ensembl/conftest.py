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
#from builtins import open
from builtins import str as text

import sys
#import os

import pytest
# import numpy as np

# from genometools.expression import ExpGene, ExpGenome, ExpMatrix

ALL = set("darwin linux cygwin win32".split())

@pytest.fixture(scope='session')
def my_download_dir(tmpdir_factory):
    download_dir = text(tmpdir_factory.mktemp('ensembl_download',
                                              numbered=False))
    return download_dir

#@pytest.fixture(scope='session'):
#def my_gene_annotation_file(my_data_pypath):
#    gene_annotation_file = text(my_gene_annotation_file.join(
#        'gene_annotation'))

def pytest_runtest_setup(item):
    # modified from: http://doc.pytest.org/en/latest/example/markers.html
    if isinstance(item, item.Function):
        plat = sys.platform
        if plat in ['linux2', 'linux3']:
            plat = 'linux'
        if not hasattr(item.obj, plat):
            if ALL.intersection(set(item.obj.__dict__)):
                pytest.skip("cannot run on platform %s" % plat)