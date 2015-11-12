#!/usr/bin/env python2.7

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

"""Command-line interface for extracting lists of protein-coding genes.

The script contained in the `main` function reads an Ensembl gene
annotation file and extracts a list of all protein-coding genes contained in
that file.

Examples
--------

Extract the protein coding genes from the human Ensembl v82 gene annotations,
downloaded from the
`Ensembl FTP server <ftp://ftp.ensembl.org/pub/release-82/gtf/homo_sapiens/>`_:

.. code-block:: bash

    $ extract_protein_coding_genes.py \\
        -a Homo_sapiens.GRCh38.82.gtf.gz \\
        -o protein_coding_genes_human.tsv

"""

import sys
import os
import gzip
import csv
import re
import gzip
import argparse
import logging
from collections import Counter

from genometools import misc

def get_argument_parser():
    """Function to obtain the argument parser.

    Returns
    -------
    A fully configured `argparse.ArgumentParser` object.

    Notes
    -----
    This function is used by the `sphinx-argparse` extension for sphinx.

    """

    parser = argparse.ArgumentParser(description='Extracts all protein-coding genes from an Ensembl gene annotation (GTF) file.')

    parser.add_argument('-a','--annotation-file',default='-',help='Path of Ensembl gene annotation file (in GTF format). Use "-" to read from stdin.')
    parser.add_argument('-o','--output-file',required=True,help='Path of output file.')
    parser.add_argument('-s','--species',choices=['human','mouse','fly','worm','fish','yeast'],default='human',help='Species for which to extract genes.')
    parser.add_argument('-c','--chromosome-pattern',required=False,default=None,
        help="""Regular expression used to determine which chromosomes to include.
                (Takes precedence over the "species" parameter when set.)""")
    parser.add_argument('-f','--field-name',default='gene',help='Only lines in the GTF file that contain this value in the third column are included.')
    parser.add_argument('-l','--log-file',default=None,help='Path of log file. If not specified, print to stdout.')
    parser.add_argument('-v','--verbose',action='store_true',help='Enable verbose output.')

    return parser

def parse_attributes(s):
    """ Parses the ``attribute`` string of a GFF/GTF annotation.

    Parameters
    ----------
    s : str
        The attribute string.

    Returns
    -------
    Dict
        A dictionary containing attribute name/value pairs.

    Notes
    -----
    The ``attribute`` string is the 9th field of each annotation (row),
    as described in the
    `GTF format specification <http://mblab.wustl.edu/GTF22.html>`_.

    """
    attr_sep = re.compile(r"(?<!\\)\s*;\s*") # use negative lookbehind to make sure we don't split on escaped semicolons ("\;")
    attr = {}
    atts = attr_sep.split(s)
    for a in atts:
        #print a
        kv = a.split(' ')
        if len(kv) == 2:
            k,v = kv
            v = v.strip('"')
            attr[k] = v
    return attr 

def main(args=None):
    """Extract protein-coding genes and store in tab-delimited text file.

    This function is the main function of the extract_protein_coding_genes.py
    script, which parses a GTF file, extract information about all protein-
    coding genes, and writes the results to a tab-delimited text file. Each
    row in the output file corresponds to one protein-coding gene. The
    columns of the output file are: 1) gene symbol, 2) chromosome name,
    3) Ensembl ID. Some genes, such as those in the
    `Pseudoautosomal region <https://en.wikipedia.org/wiki/Pseudoautosomal_region>`_,
    have annotations for multiple chromosomes, and/or are associated with
    multiple Ensembl IDs. In those cases, columns 2) and/or 3) contain
    all values, separated by a comma.

    Parameters
    ----------
    args: argparse.Namespace object, optional
        The argument values. If not specified, the values will be obtained by
        parsing the command line arguments using the `argparse` module.

    Returns
    -------
    int
        Exit code (0 if no error occurred).
 
    """

    chromosome_patterns = {\
            'human': r'(?:\d\d?|MT|X|Y)$',\
            'mouse': r'(?:\d\d?|MT|X|Y)$',\
            'fly': r'(?:2L|2R|3L|3R|4|X|Y|dmel_mitochondrion_genome)$',\
            'worm': r'(?:I|II|III|IV|V|X|MtDNA)$',\
            'fish': r'(?:\d\d?|MT)$',\
            'yeast': r'(?:I|II|III|IV|V|VI|VII|VIII|IX|X|XI|XII|XIII|XIV|XV|XVI|Mito)$'}

    if args is None:
        parser = get_argument_parser()
        args = parser.parse_args()

    input_file = args.annotation_file
    species = args.species
    chrom_pat = args.chromosome_pattern
    field_name = args.field_name
    log_file = args.log_file
    verbose = args.verbose

    # configure logger
    log_level = logging.INFO
    if verbose:
        log_level = logging.DEBUG

    log_format = '[%(asctime)s] %(levelname)s: %(message)s'
    log_datefmt = '%Y-%m-%d %H:%M:%S'
    # when filename is not None, "stream" parameter is ignored (see https://docs.python.org/2/library/logging.html#logging.basicConfig)
    logging.basicConfig(filename=log_file,stream=sys.stdout,level=log_level,format=log_format,datefmt=log_datefmt)
    logger = logging.getLogger()

    if chrom_pat is None:
        chrom_pat = re.compile(chromosome_patterns[species])
    else:
        chrom_pat = re.compile(chrom_pat)

    #if exclude_chromosomes:
    #   print "Excluding chromosomes %s..." %(', '.join(sorted(exclude_chromosomes)))
    #   sys.stdout.flush()
    logger.info('Regular expression used for filtering chromosome names: %s',chrom_pat.pattern)

    # for statistics
    types = Counter()
    sources = Counter()

    # primary information
    genes = Counter()
    gene_chroms = dict()
    gene_ids = dict()

    # secondary information
    genes2 = Counter()
    polymorphic = set()

    # list of chromosomes
    chromosomes = set()

    i = 0
    missing = 0
    excluded_chromosomes = set()
    logger.info('Parsing data...')
    with misc.open_plain_or_gzip(input_file) if input_file != '-' else sys.stdin as fh:
        #if i >= 500000: break
        reader = csv.reader(fh,dialect='excel-tab')
        for l in reader:
            i += 1
            #if i % int(1e5) == 0:
            #   print '\r%d...' %(i), ; sys.stdout.flush() # report progress

            if len(l) > 1 and l[2] == field_name:

                # Note: Older Ensembl GTF files sometimes have a leading space in field #9
                attr = parse_attributes(l[8].lstrip(' '))

                type_ = attr['gene_biotype']
                if type_ not in ['protein_coding','polymorphic_pseudogene']:
                    continue

                chrom = l[0]

                # test whether chromosome is valid
                m = chrom_pat.match(chrom)
                if m is None:
                    excluded_chromosomes.add(chrom)
                    continue

                chromosomes.add(m.group())

                source = l[1]
                id_ = attr['gene_id']
                try:
                    name = attr['gene_name']
                except KeyError as e:
                    missing += 1
                    continue

                # store gene name
                genes[name] += 1

                # store Ensemble ID
                try:
                    gene_ids[name].add(id_)
                except KeyError:
                    gene_ids[name] = set([id_])

                # store chromsome
                try:
                    gene_chroms[name].add(chrom)
                except KeyError:
                    gene_chroms[name] = set([chrom])

                # record some statistics
                sources[source] += 1
                types[type_] += 1
                if type_ == 'polymorphic_pseudogene':
                    polymorphic.add(name)
                    genes2[name] += 1

    logger.info('done (parsed %d lines).',i)

    logger.info('')
    logger.info('Gene chromosomes (%d):',len(chromosomes))
    logger.info('\t' + ', '.join(sorted(chromosomes)))
    logger.info('')
    logger.info('Excluded chromosomes (%d):',len(excluded_chromosomes))
    logger.info('\t' + ', '.join(sorted(excluded_chromosomes)))

    if missing > 0:
        logger.info('')
        logger.info('# Genes without names: %d',missing)

    def print_counter(C):
        for k in sorted(C.keys(),key=lambda x:-C[x]):
            logger.info('\t%s: %d' %(k,C[k]))

    logger.info('')
    logger.info('Gene sources:')
    print_counter(sources)

    logger.info('')
    logger.info('Gene types:')
    print_counter(types)

    redundant_genes = sorted(g for g in genes if genes[g]>1)
    logger.info('')
    logger.info('# Genes with redundant annotations: %d',len(redundant_genes))

    logger.info('')
    logger.info('Polymorphic pseudogenes (%d): %s',len(polymorphic), ', '.join('%s (%d)' %(g,genes2[g]) for g in sorted(polymorphic)))

    logger.info('')
    logger.info('Total protein-coding genes: %d', len(genes))

    with open(args.output_file,'w') as ofh:
        writer = csv.writer(ofh,dialect='excel-tab',lineterminator='\n',quoting=csv.QUOTE_NONE)
        for name in sorted(genes):
            chroms = ','.join(sorted(gene_chroms[name]))
            ids = ','.join(sorted(gene_ids[name]))
            writer.writerow([name,chroms,ids])

    return 0

if __name__ == '__main__':
    return_code = main()
    sys.exit(return_code)
