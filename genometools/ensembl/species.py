# Copyright (c) 2015, 2016 Florian Wagner
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
_oldstr = str
from builtins import *

from collections import OrderedDict


SPECIES_NAMES = OrderedDict([
    ('Homo_sapiens', 'human'),
    ('Mus_musculus', 'mouse'),
    ('Danio_rerio', 'zebrafish'),
    ('Drosophila_melanogaster', 'fly'),
    ('Caenorhabditis_elegans', 'worm'),
    ('Saccharomyces_cerevisiae', 'yeast'),
])
"""Dictionary mapping scientific species names to their common names."""

SPECIES_SCIENTIFIC = OrderedDict(
    [b, a] for a, b in SPECIES_NAMES.items()
)
"""Reverse mapping of ``species_names``."""

SPECIES_CHROMPAT = OrderedDict([
    ('human', r'(?:\d\d?|MT|X|Y)$'),
    ('mouse', r'(?:\d\d?|MT|X|Y)$'),
    ('fly', r'(?:2L|2R|3L|3R|4|X|Y|dmel_mitochondrion_genome)$'),
    ('worm', r'(?:I|II|III|IV|V|X|MtDNA)$'),
    ('zebrafish', r'(?:\d\d?|MT)$'),
    ('yeast', r'(?:I|II|III|IV|V|VI|VII|VIII|IX|X|XI|XII|XIII|XIV|XV|XVI|'
              r'Mito)$')
])
"""Regular expressions implicitly defining "valid" Ensembl chromosome names for
each species."""
