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

"""Tests for the `ensembl.annotation` module."""

import hashlib

from pandas.util import hash_pandas_object

from genometools import ensembl


def test_get_genes(my_gene_annotation_file):
    """Test the extraction of genes from gene annotation (GTF) files."""

    # test protein-coding genes
    genes = ensembl.get_protein_coding_genes(my_gene_annotation_file)
    assert len(genes) == 4
    h = hashlib.md5(hash_pandas_object(genes).values.tobytes()).hexdigest()
    assert h == 'a72941001a18af44ffdd555c4fe30bd8'
