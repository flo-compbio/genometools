from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import *

from collections import OrderedDict

species_names = OrderedDict([
        ('Homo_sapiens', 'human'),
        ('Mus_musculus', 'mouse'),
        ('Danio_rerio', 'zebrafish'),
        ('Drosophila_melanogaster', 'fly'),
        ('Caenorhabditis_elegans', 'worm'),
        ('Saccharomyces_cerevisiae', 'yeast'),
])
"""Dictionary mapping scientific species names to their common names."""

species_scientific = OrderedDict(
        [b, a] for a, b in species_names.items()
)
"""Reverse mapping of ``species_names``."""

species_chrompat = OrderedDict([
        ('human', r'(?:\d\d?|MT|X|Y)$'),
        ('mouse', r'(?:\d\d?|MT|X|Y)$'),
        ('fly', r'(?:2L|2R|3L|3R|4|X|Y|dmel_mitochondrion_genome)$'),
        ('worm', r'(?:I|II|III|IV|V|X|MtDNA)$'),
        ('zebrafish', r'(?:\d\d?|MT)$'),
        ('yeast', r'(?:I|II|III|IV|V|VI|VII|VIII|IX|X|XI|XII|XIII|XIV|XV|XVI|Mito)$')
])
"""Regular expressions implicitly defining "valid" Ensembl chromosome names for
each species."""
