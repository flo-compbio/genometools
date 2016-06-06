# Copyright (c) 2016 Florian Wagner
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

"""Heat map annotation classes.

"""

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import *

import logging

logger = logging.getLogger(__name__)


class HeatMapAnnotation(object):

    # TODO: docstrings, __str__, __repr__, ...

    def __init__(self, key, color, label=None):

        # key can be anything allowed as a key in a pandas index.
        assert isinstance(color, str)
        assert isinstance(label, str)

        self.key = key
        self.color = color
        self.label = label


class HeatMapGeneAnnotation(HeatMapAnnotation):

    # TODO: docstrings, __str__, __repr__, ...

    def __init__(self, gene, color, label=None):
        HeatMapAnnotation.__init__(self, gene, color, label)

    @property
    def gene(self):
        """Alias for `HeatMapAnnotation.key`."""
        return self.key

    @gene.setter
    def gene(self, value):
        self.key = value


class HeatMapSampleAnnotation(HeatMapAnnotation):

    # TODO: docstrings, __str__, __repr__, ...

    def __init__(self, sample, color, label=None):
        HeatMapAnnotation.__init__(self, sample, color, label)

    @property
    def sample(self):
        """Alias for `HeatMapAnnotation.key`."""
        return self.key

    @sample.setter
    def sample(self, value):
        self.key = value