# Copyright (c) 2016 Florian Wagner
#
# This file is part of GO-PCA.
#
# GO-PCA is free software: you can redistribute it and/or modify
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

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
# from builtins import *
from builtins import open
from builtins import str as text

# import os
import shutil

import pytest
import requests
import logging

logging.basicConfig(level=logging.INFO)

logger = logging.getLogger(__name__)

from genometools import misc
from genometools.ontology import GOTerm, GeneOntology
from genometools.expression import ExpGenome

def download_file(url, path):
    r = requests.get(url, stream=True)
    if r.status_code == 200:
        with open(path, 'wb') as f:
            r.raw.decode_content = True
            shutil.copyfileobj(r.raw, f)


@pytest.fixture(scope='session')
def my_data_pypath(tmpdir_factory):
    pypath = tmpdir_factory.mktemp('gopca_data', numbered=False)
    return pypath


@pytest.fixture(scope='session')
def my_genome_file(my_data_pypath):
    logger.info('Starting download of genome file...')
    url = r'https://www.dropbox.com/s/vvhz0jt2ly9hzsz/protein_coding_genes_human_ensembl83.tsv?dl=1'
    path = text(my_data_pypath.join('protein_coding_genes_human_ensembl83.tsv'))
    download_file(url, path)
    return path


@pytest.fixture(scope='session')
def my_genome(my_genome_file):
    genome = ExpGenome.read_tsv(my_genome_file)
    return genome


@pytest.fixture(scope='session')
def my_gene_ontology_file(my_data_pypath):
    logger.info('Starting download of gene ontology file...')
    url = r'https://www.dropbox.com/s/rvox0wi96162it2/go-basic_2016-01-18.obo?dl=1'
    path = text(my_data_pypath.join('go-basic_2016-01-18.obo'))
    download_file(url, path)
    return path


@pytest.fixture(scope='session')
def my_gene_ontology(my_gene_ontology_file):
    gene_ontology = GeneOntology.read_obo(my_gene_ontology_file)
    return gene_ontology


@pytest.fixture(scope='session')
def my_goa_file(my_data_pypath):
    logger.info('Starting download of GO annotation file...')
    url = r'ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/HUMAN/goa_human.gaf.gz'
    path = text(my_data_pypath.join('goa_human.gaf.gz'))
    misc.ftp_download(url, path)
    return path


@pytest.fixture(scope='session')
def my_go_term():
    go_term = GOTerm('GO:0000000', 'regulation of test process',
                  'biological_process', 'This is a test GO term.')
    return go_term