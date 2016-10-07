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
from builtins import open
from builtins import str as text

import pytest


@pytest.fixture(scope='session')
def my_temp_dir(tmpdir_factory):
    temp_dir = tmpdir_factory.mktemp('gopca_output', numbered=False)
    return temp_dir


@pytest.fixture(scope='session')
def my_checksum_file(my_temp_dir):
    checksum_file = text(my_temp_dir.join('checksum_file.txt'))
    with open(checksum_file, 'w', encoding='UTF-8') as ofh:
        ofh.write('Hello world!')
    print('Test:', open(checksum_file).read())
    return checksum_file