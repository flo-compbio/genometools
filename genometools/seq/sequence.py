# Copyright (c) 2015 Florian Wagner
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

"""Module containing the `Sequence` class."""

class Sequence(object):
    """A nucleotide sequence."""

    def __init__(self, name, seq):
        assert isinstance(name, str)
        assert isinstance(seq, str)
        self.name = name
        self.seq = seq

    def __repr__(self):
        return '<%s "%s" (length=%d; seq_hash=%d)>' %(self.__class__.__name__,
                self.name, self.length, self.seq_hash)

    def __str__(self):
        return '<%s "%s" (length=%d)>' \
                %(self.__class__.__name__, self.name, self.length)

    def __eq__(self, other):
        if self is other:
            return True
        elif type(self) != type(other):
            return False
        else:
            return repr(self) == repr(other)

    def __hash__(self):
        data = []
        data.append(self.name)
        data.append(self.seq)
        return hash(tuple(data))

    @property
    def length(self):
        return len(self.seq)

    @property
    def seq_hash(self):
        return hash(self.seq)

    def append_fasta(self, ofh, linewidth = 70):
        ofh.write('>%s\n' %(self.name))
        start = 0
        while start < self.length - linewidth:
            ofh.write('%s\n' %(self.seq[start:(start+linewidth)]))
            start += linewidth
        ofh.write('%s\n' %(self.seq[start:]))
