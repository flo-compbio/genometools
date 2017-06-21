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

"""Functions for working with Ensembl cDNA data."""

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
_oldstr = str
from builtins import *

import os
import gzip
from collections import OrderedDict
import logging

import pandas as pd
from Bio import SeqIO


_LOGGER = logging.getLogger(__name__)


def _transform_chrom(chrom):
    """Helper function to obtain specific sort order."""
    try:
        c = int(chrom)
    except:
        if chrom in ['X', 'Y']:
            return chrom
        elif chrom == 'MT':
            return '_MT'  # sort to the end
        else:
            return '__' + chrom # sort to the very end
    else:
        # make sure numbered chromosomes are sorted numerically
        return '%02d' % c


def get_chromosome_lengths(fasta_file, fancy_sort=True):
    """Extract chromosome lengths from genome FASTA file."""
    chromlen = []
    with gzip.open(fasta_file, 'rt', encoding='ascii') as fh:
        fasta = SeqIO.parse(fh, 'fasta')
        for i, f in enumerate(fasta):
            chromlen.append((f.id, len(f.seq)))
            _LOGGER.info('Processed chromosome "%s"...', f.id)
            #print(dir(f))
            #if i == 1: break
            
    # convert to pandas Series
    chromlen = pd.Series(OrderedDict(chromlen))
    chromlen.index.name = 'Chromosome'
    chromlen.name = 'Length'

    if fancy_sort:
        # sort using fancy ordering
        chrom_for_sorting = chromlen.index.to_series().apply(_transform_chrom)
        a = chrom_for_sorting.argsort(kind='mergesort')
        chromlen = chromlen.iloc[a]

    return chromlen
