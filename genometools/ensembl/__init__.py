species_chrompat = {\
        'human': r'(?:\d\d?|MT|X|Y)$',\
        'mouse': r'(?:\d\d?|MT|X|Y)$',\
        'fly': r'(?:2L|2R|3L|3R|4|X|Y|dmel_mitochondrion_genome)$',\
        'worm': r'(?:I|II|III|IV|V|X|MtDNA)$',\
        'fish': r'(?:\d\d?|MT)$',\
        'yeast': r'(?:I|II|III|IV|V|VI|VII|VIII|IX|X|XI|XII|XIII|XIV|XV|XVI|Mito)$'}
"""Regular expressions implicitly defining "valid" chromosome names for each
species."""
