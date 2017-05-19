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

#from __future__ import (absolute_import, division,
#                        print_function, unicode_literals)
#from builtins import *

import os
import pkg_resources

__version__ = pkg_resources.require('genometools')[0].version

_root = os.path.abspath(os.path.dirname(__file__))

from . import basic
from . import ebi
from . import enrichment
from . import ensembl
from . import expression
from . import gdc
from . import misc
from . import ncbi
from . import ontology
from . import rnaseq
from . import seq

# __all__ = ['misc', 'gtf', 'ensembl', 'rnaseq', 'rnaseq']
