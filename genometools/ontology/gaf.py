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

"""Module containing the `GOACollection` class."""

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
_oldstr = str
from builtins import *

import sys
import re
import hashlib
import logging
from collections import OrderedDict, Iterable

import pandas as pd
import numpy as np

from . import GOTerm, GeneOntology, GOAnnotation
from .. import misc
from ..basic import GeneSet, GeneSetCollection

logger = logging.getLogger(__name__)


def parse_gaf(path_or_buffer, gene_ontology, genome=None,
              db=None, ev_codes=None):
    """Parse a GAF 2.1 file containing GO annotations.
    
    Parameters
    ----------
    path_or_buffer : str or buffer
        The GAF file.
    gene_ontology : `GeneOntology`
        The Gene Ontology.
    genome : `expression.ExpGenome`
        The genome.
    db : str, optional
        Select only annotations with this "DB"" value.
    ev_codes : str or set of str, optional
        Select only annotations with this/these evidence codes.
    
    Returns
    -------
    list of `GOAnnotation`
        The list of GO annotations.
    """
    #if path == '-':
    #    path = sys.stdin

    assert isinstance(gene_ontology, GeneOntology)
    if db is not None:
        assert isinstance(db, (str, _oldstr))
    if (ev_codes is not None) and ev_codes:
        assert isinstance(ev_codes, (str, _oldstr)) or \
                isinstance(ev_codes, Iterable)

    if isinstance(ev_codes, str):
        ev_codes = set([ev_codes])
    elif (ev_codes is not None) and ev_codes:
        ev_codes = set(ev_codes)
    else:
        ev_codes = None

    # open file, if necessary
    if isinstance(path_or_buffer, (str, _oldstr)):
        buffer = misc.gzip_open_text(path_or_buffer, encoding='ascii')
    else:
        buffer = path_or_buffer

    # use pandas to parse the file quickly
    df = pd.read_csv(buffer, sep='\t', comment='!', header=None, dtype=_oldstr)

    # replace pandas' NaNs with empty strings
    df.fillna('', inplace=True)

    # exclude annotations with unknown Gene Ontology terms
    all_go_term_ids = set(gene_ontology._term_dict.keys())
    sel = df.iloc[:, 4].isin(all_go_term_ids)
    logger.info(
        'Ignoring %d / %d annotations (%.1f %%) with unknown GO terms.',
        (~sel).sum(), sel.size, 100*((~sel).sum()/float(sel.size)))
    df = df.loc[sel]

    # filter rows for genome
    if genome is not None:
        all_genes = set(genome.gene_names)
        sel = df.iloc[:, 2].isin(all_genes)
        logger.info(
            'Ignoring %d / %d annotations (%.1f %%) with unknown genes.',
            (~sel).sum(), sel.size, 100*((~sel).sum()/float(sel.size)))
        df = df.loc[sel]

    # filter rows for DB value
    if db is not None:
        sel = (df.iloc[:, 0] == db)
        logger.info(
            'Excluding %d / %d annotations (%.1f %%) with wrong DB values.',
            (~sel).sum(), sel.size, 100*((~sel).sum()/float(sel.size)))
        df = df.loc[sel]

    # filter rows for evidence value
    if ev_codes is not None:
        sel = (df.iloc[:, 6].isin(ev_codes))
        logger.info(
            'Excluding %d / %d annotations (%.1f %%) based on evidence code.',
            (~sel).sum(), sel.size, 100*((~sel).sum()/float(sel.size)))
        df = df.loc[sel]

    # convert each row into a GOAnnotation object
    go_annotations = []
    for i, l in df.iterrows():
        ann = GOAnnotation.from_list(gene_ontology, l.tolist())
        go_annotations.append(ann)
    logger.info('Read %d GO annotations.', len(go_annotations))

    return go_annotations


def get_goa_gene_sets(go_annotations):
    """Generate a list of gene sets from a collection of GO annotations.

    Each gene set corresponds to all genes annotated with a certain GO term.
    """
    go_term_genes = OrderedDict()
    term_ids = {}
    for ann in go_annotations:
        term_ids[ann.go_term.id] = ann.go_term
        try:
            go_term_genes[ann.go_term.id].append(ann.db_symbol)
        except KeyError:
            go_term_genes[ann.go_term.id] = [ann.db_symbol]
    
    go_term_genes = OrderedDict(sorted(go_term_genes.items()))
    gene_sets = []
    for tid, genes in go_term_genes.items():
        go_term = term_ids[tid]
        gs = GeneSet(id=tid, name=go_term.name, genes=genes,
                     source='GO',
                     collection=go_term.domain_short,
                     description=go_term.definition)
        gene_sets.append(gs)
    gene_sets = GeneSetCollection(gene_sets)
    return gene_sets