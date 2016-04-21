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

"""Module containing the `FastaReader` class."""

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import *

from . import Sequence


class FastaReader(object):
    """A FASTA file reader."""

    def __init__(self, fh):
        self.fh = fh
        self.cur_name = None

    def __iter__(self):
        return self

    def __next__(self):
        started = False
        name = None
        if self.cur_name is not None:
            started = True
            name = self.cur_name

        seq = []

        while True:

            try:
                l = self.fh.next()
            except StopIteration:
                if started:
                    # we've read a sequence and now hit the end of the file
                    self.cur_name = None
                    break
                else:
                    # we haven't read a sequence, so the iterator is done
                    raise StopIteration

            if l[0] == '>':
                if started:
                    # this is already the beginning of the next sequence,
                    # => store the name of the next sequence, so that the
                    #    information from this line isn't lost
                    self.cur_name = l[1:-1]
                    break

                # this is the beginning of the sequence
                started = True
                name = l[1:-1]

            elif started:
                if l[0] == '\n':  # we've hit an empty line (= sequence end)
                    self.cur_name = None
                    break
                seq.append(l[:-1])
            
        return Sequence(name, ''.join(seq))
