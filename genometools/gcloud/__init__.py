# Copyright (c) 2017 Florian Wagner
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

#from __future__ import (absolute_import, division,
#                        print_function, unicode_literals)
#from builtins import *

#import os
#import pkg_resources

# __version__ = pkg_resources.require('genometools')[0].version

#_root = os.path.abspath(os.path.dirname(__file__))

from .service_account import *
from .operations import *

from . import storage
from . import compute
from . import tasks
