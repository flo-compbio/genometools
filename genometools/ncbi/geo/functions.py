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

"""Functions for files from the NCBI GEO database."""

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import *

import os
from collections import Iterable

import unicodecsv as csv
import numpy as np

from ... import misc


def get_series_matrix_url(acc):
    assert isinstance(acc, str)
    return 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/matrix/' \
           '%s_series_matrix.txt.gz' \
            % (acc[:5] + 'nnn', acc, acc)


def get_raw_data_url(acc):
    assert isinstance(acc, str)
    return 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/suppl/%s_RAW.tar' \
            % (acc[:5] + 'nnn', acc, acc)


def read_series_matrix(path):
    """Extracts sample information from a GEO series matrix file.

    The extracted information consists of sample accessions, names, and raw
    data URLs.

    Parameters
    ----------
    path: str
        The series matrix file. Can be gzip'ed.

    Returns
    -------
    3-tuple (accessions, titles, raw data URLs)
    """
    assert isinstance(path, str)

    accessions = None
    names = None
    celfile_urls = None
    with misc.smart_open_read(path, mode='rb', try_gzip=True) as fh:
        reader = csv.reader(fh, dialect='excel-tab')
        for l in reader:
            if not l:
                continue
            if l[0] == '!Sample_geo_accession':
                accessions = l[1:]
            elif l[0] == '!Sample_title':
                names = l[1:]
            elif l[0] == '!Sample_supplementary_file':
                celfile_urls = l[1:]
            elif l[0] == '!series_matrix_table_begin':
                # we've read the end of the section containing metadata
                break

    return accessions, names, celfile_urls


def write_sample_sheet(output_file, accessions, names, celfile_urls, sel=None):
    """Generate a sample sheet in tab-separated text format.

    The columns contain the following sample attributes:
    1) accession
    2) name
    3) CEL file name
    4) CEL file URL

    Parameters
    ----------
    output_file: str
        The path of the output file.
    accessions: list or tuple of str
        The sample accessions.
    names: list or tuple of str
        The sample names.
    celfile_urls: list or tuple of str
        The sample CEL file URLs.
    sel: Iterable, optional
        A list of sample indices to include. If None, all samples are included.
        [None]

    Returns
    -------
    None
    """
    assert isinstance(output_file, str)
    assert isinstance(accessions, (list, tuple))
    for acc in accessions:
        assert isinstance(acc, str)
    assert isinstance(names, (list, tuple))
    for n in names:
        assert isinstance(n, str)
    assert isinstance(celfile_urls, (list, tuple))
    for u in celfile_urls:
        assert isinstance(u, str)
    if sel is not None:
        assert isinstance(sel, Iterable)
        for i in sel:
            assert isinstance(i, (int, np.integer))

    with open(output_file, 'wb') as ofh:
        writer = csv.writer(ofh, dialect='excel-tab',
                            lineterminator=os.linesep,
                            quoting=csv.QUOTE_NONE)
        # write header
        writer.writerow(['Accession', 'Name', 'CEL file name', 'CEL file URL'])
        n = len(list(names))
        if sel is None:
            sel = range(n)
        for i in sel:
            cf = celfile_urls[i].split('/')[-1]
            # row = [accessions[i], names[i], cf, celfile_urls[i]]
            writer.writerow([accessions[i], names[i], cf, celfile_urls[i]])
