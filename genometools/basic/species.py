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

"""Module containing the `Species` class.

"""

import logging

logger = logging.getLogger(__name__)

class Species(object):
    """A species.

    Parameters
    ----------
    name: str
        See :attr:`name` attribute.
    common_name: str
        See :attr:`short_name` attribute.

    Attributes
    ----------
    name: str
        The scientific name of the species (e.g., "Homo sapiens").
    common_name: str, optional
        The common name of the species (e.g., "human")
    """
    def __init__(self, name, common_name = None):
        self.name = name
        self.common_name = common_name

    def __repr__(self):
        return '<Species "%s">' %(self.name)

    def __str__(self):
        common = ''
        if self.common_name is not None:
            common = ' (%s)' %(self.common_name)
        return '<Species "%s"%s>' %(self.name, common)

    def __hash__(self):
        data = []
        data.append(self.name)
        return hash(frozenset(data))

    def __eq__(self,other):
        if type(self) != type(other):
            return False
        else:
            return repr(self) == repr(other)
