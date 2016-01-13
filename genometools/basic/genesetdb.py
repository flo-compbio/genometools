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

"""Module containing the `GeneSetDB` class.

"""

import os
import logging
from collections import Iterable, OrderedDict
#import codecs

import xmltodict
import unicodecsv as csv

from genometools.basic import GeneSet

logger = logging.getLogger(__name__)

class GeneSetDB(object):
    """A gene set database."""

    def __init__(self, gene_sets):
        
        assert isinstance(gene_sets, Iterable)
        for gs in gene_sets:
            assert isinstance(gs, GeneSet)

        self.gene_sets = OrderedDict([gs.id, gs] for gs in gene_sets)

    def __repr__(self):
        return '<%s (n=%d; hash=%d)>' \
                %(self.__class__, self.n, hash(self))

    def __hash__(self):
        return(hash(frozenset(self.gene_sets.values())))

    @property
    def n(self):
        return len(self.gene_sets)

    @classmethod
    def read_tsv(cls, path):
        """Read the database from a tab-delimited text file."""
        gene_sets = []
        with open(path, 'rb') as fh:
            reader = csv.reader(fh, dialect = 'excel-tab')
            for l in reader:
                gs = GeneSet.from_list(l)
                gene_sets.append(gs)
        return cls(gene_sets)

    def write_tsv(self, path):
        """Write the database to a tab-delimited text file."""
        with open(path, 'wb') as ofh:
            writer = csv.writer(ofh, dialect = 'excel-tab',
                quoting = csv.QUOTE_NONE, lineterminator = os.linesep)
            for gs in self.gene_sets.itervalues():
                writer.writerow(gs.to_list())

    @classmethod
    def read_msigdb_xml(cls, path, entrez2gene, species = None):
        """Read the complete MSigDB database from an XML file.

        The XML file can be downloaded from here:
        http://software.broadinstitute.org/gsea/msigdb/download_file.jsp?filePath=/resources/msigdb/5.0/msigdb_v5.0.xml
        """

        logger.debug('Path: %s', path)
        logger.debug('entrez2gene type: %s', str(type(entrez2gene)))

        i = [0]
        gene_sets = []

        total_gs = [0]
        total_genes = [0]

        species_excl = [0]
        unknown_entrezid = [0]

        src = 'MSigDB'

        def handle_item(pth, item):

            total_gs[0] += 1
            data = pth[1][1]

            spec = data['ORGANISM']
            # filter by species
            if species is not None and spec != species:
                species_excl[0] += 1
                return True

            id_ = data['SYSTEMATIC_NAME']
            name = data['STANDARD_NAME']
            coll = data['CATEGORY_CODE']
            desc = data['DESCRIPTION_BRIEF']
            entrez = data['MEMBERS_EZID'].split(',')

            genes = []
            for e in entrez:
                total_genes[0] += 1
                try:
                    genes.append(entrez2gene[e])
                except KeyError:
                    unknown_entrezid[0] += 1

            if not genes:
                logger.warning('Gene set "%s" (%s) has no known genes!',
                        name, id_)
                return True

            gs = GeneSet(name, genes, id_ = id_, source = src,
                    collection = coll, description = desc)
            gene_sets.append(gs)
            i[0] += 1
            return True

        with open(path, 'rb') as fh:
            print type(fh)
            xmltodict.parse(fh.read(), encoding = 'UTF-8', item_depth = 2,
                    item_callback = handle_item)

        if species_excl[0] > 0:
            kept = total_gs[0] - species_excl[0]
            perc = 100 * (kept / float(total_gs[0]))
            logger.info('%d of all %d gene sets (%.1f %%) belonged to the ' +
                    'specified species.', kept, total_gs[0], perc)

        if unknown_entrezid[0] > 0:
            unkn = unknown_entrezid[0]
            #known = total_genes[0] - unknown_entrezid[0]
            perc = 100 * (unkn / float(total_genes[0]))
            logger.warning('%d of a total of %d genes (%.1f %%) had an ' +
                    'unknown Entrez ID.', unkn, total_genes[0], perc)

        logger.info('Parsed %d entries, resulting in %d gene sets.',
                total_gs[0], len(gene_sets))

        return cls(gene_sets)
