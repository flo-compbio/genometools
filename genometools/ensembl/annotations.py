# Copyright (c) 2016 Florian Wagner
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

import os
import time
import re
import logging

import pandas as pd
    
from genometools import misc
from genometools import gtf

logger = logging.getLogger(__name__)


def get_protein_coding_genes(
        path_or_buffer, chunksize=100000,
        chromosome_pattern=r'(?:\d\d?|MT|X|Y)$',
        include_polymorphic_pseudogenes=True,
        only_manual=False,
        remove_duplicates=True,
        fancy_sorting=True):
    r"""Get list of all protein-coding genes based on Ensembl GTF file.
    
    Parameters
    ----------
    path_or_buffer : str or buffer
        The GTF file (either the file path or a buffer)
    chromosome_pattern : str, optional
        Regular expression specifying valid chromosomes. [r'(?:\d\d?|MT|X|Y)$']
    include_polymorphic_pseudogene : bool, optional
        Whether to include genes annotated as "polymorphic pseudogenes"?
    only_manual : bool, optional
        Whether to exclude annotations with source "ensembl", which
        are based only on an automatic annotation pipeline. [True]
    remove_duplicates : bool, optional
        Whether to remove duplicate annotations, i.e. those with different
        Ensembl IDs for the same gene. [True]
    fancy_sorting : bool, optional
        Whether to sort chromosomes numerically, with "X", "Y", and "MT" at the
        end. 

    Returns
    -------
    `pandas.DataFrame`
        Table with rows corresponding to protein-coding genes.

    Notes
    -----
    
    Annotation sources and redundant gene annotations
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    According to the Ensembl website (1), the Ensembl gene annotation
    GTF files for human, mouse, zebrafish, rat and pig essentially
    contain two sets of annotations:
    
    One set consists of all annotations with the "ensembl"
    source annotation (column 2). These annotations are the product of
    the automated Ensembl "genebuild" pipeline.
    
    The other set consists of genes that are manually annotated by
    the HAVANA team (source "havana"), some of which have been merged with the
    automatic annotations (source "ensembl_havana").
    
    There seems to be no overlap between genes annotated with "havana" and
    "ensembl_havana" sources, respectively. However, there are a few genes for
    which only annotations with source "ensembl" exist.
    
    Our policy is therefore to prefer annotations with source "ensembl_havana"
    and "havana" over those with source "ensembl", and to only keep annotations
    with source "ensembl" if there are no manually curated alternative
    annotations.
    
    A special case is represented by mitochondrial genes, which always have the
    source "insdc".
    
    (1) see http://www.ensembl.org/Help/Faq?id=152


    Removal of duplicates
    ~~~~~~~~~~~~~~~~~~~~~
    
    Unfortunately, the Ensembl gene annotations contain duplicates for a
    handful of genes. For example, for MATR3, there are ENSG00000015479 and
    ENSG00000280987, both of type
    "ensembl_havana". There seems to be no clear criterion by which we could
    rationally and automatically choose one ID over the other, at least based
    on information contained
    in the GTF file.
    
    We therefore remove duplicates according to following policy:
    - For genes on '+' strand, keep the gene with the left-most starting
      position.
    - For genes on '-' strand, keep the gene with the right-most starting
      position.
    (In case the starting positions are equal, we keep the one that occurs
    first in the GTF file.)
    
    We would like to use the pandas.DataFrame.drop_duplicates() function for
    this. So we're temporarily reordering genes using their signed position,
    and then we're using the original index (position) to restore the original
    order.
    """
    chrompat = re.compile(chromosome_pattern)
    
    c = 0
    num_lines = 0
    num_chunks = 0
    
    t0 = time.time()
    reader = pd.read_csv(path_or_buffer, encoding='ascii', sep='\t',
                         header=None, comment='#', dtype={0: str},
                         chunksize=chunksize)
    data = []
    header = ['Gene', 'Ensembl_ID',
              'Chromosome', 'Position', 'Length',
              'Source', 'Type']
    
    valid_biotypes = set(['protein_coding'])
    if include_polymorphic_pseudogenes:
        valid_biotypes.add('polymorphic_pseudogene')
        
    valid_sources = set(['ensembl_havana', 'havana', 'insdc'])
    if not only_manual:
        valid_sources.add('ensembl')
        
    excluded_chromosomes = set()
    
    for j, df in enumerate(reader):
        num_chunks += 1
        num_lines += (df.shape[0])
        # "insdc" is required to catch the mitochondrial protein-coding genes
        sel = (df.iloc[:, 2] == 'gene') & df.iloc[:, 1].isin(valid_sources)
        # c += sel.sum()
        for i, row in df.loc[sel].iterrows():
            attr = gtf.parse_attributes(row[8].lstrip(' '))

            biotype = attr['gene_biotype']
            if biotype not in valid_biotypes:
                continue

            chrom = str(row[0])
            source = row[1]
            match = chrompat.match(chrom)
            if match is None:
                excluded_chromosomes.add(chrom)
                continue

            c += 1

            gene_name = attr['gene_name']
            ensembl_id = attr['gene_id']

            assert row[6] in ['+', '-']
            if row[6] == '+':
                pos = int(row[3])-1
            elif row[6] == '-':
                pos = -int(row[4])
            else:
                raise ValueError('Invalid strand information: %s'
                                 % str(row[6]))
            length = abs(int(row[4]) - int(row[3]))

            data.append([gene_name, ensembl_id, chrom, pos, length,
                         source, biotype])
            
    t1 = time.time()
    
    df = pd.DataFrame(columns=header, data=data)
    
    if not only_manual:
        # keep only annotations with source "ensembl"
        # if no manual annotations are available
        sel = df['Source'] == 'ensembl'
        redundant_ensembl_genes = set(df.loc[sel, 'Gene'].values) & \
                set(df.loc[~sel, 'Gene'].values)
        sel = sel & df['Gene'].isin(redundant_ensembl_genes)
        num_genes_before = df.shape[0]
        df = df.loc[~sel]
        num_genes_after = df.shape[0]
        logger.info('Removed %d gene annotations with source "ensembl" that '
                    'also had manual annotations.',
                    num_genes_before-num_genes_after)
    if remove_duplicates:
        # remove duplicate annotations (two or more Ensembl IDs for the same
        # gene)
        num_genes_before = df.shape[0]

        # sort by signed position value
        df.sort_values('Position', kind='mergesort', inplace=True)

        # remove duplicates by keeping the first occurrence
        df.drop_duplicates(['Chromosome', 'Gene'], inplace=True)

        # restore original order using the numeric index
        df.sort_index(inplace=True)

        num_genes_after = df.shape[0]
        logger.info('Removed %d duplicate gene entries',
                    num_genes_before-num_genes_after)

    # sort normally (first by chromsome, then by absolute position)
    df_sort = pd.concat([df['Chromosome'], df['Position'].abs()], axis=1)
    df_sort = df_sort.sort_values(['Chromosome', 'Position'], kind='mergesort')
    df = df.loc[df_sort.index]

    if fancy_sorting:
        # Perform "fancy sorting" of genes. Chromosomes with numbers (1-22)
        # are ordered numerically, and followed by the X, Y, and MT
        # chromosomes.
        def transform_chrom(chrom):
            try:
                c = int(chrom)
            except:
                if chrom == 'MT':
                    return '_MT'
                else:
                    return chrom
            else:
                return '%02d' % c

        chrom_for_sorting = df['Chromosome'].apply(transform_chrom)
        a = chrom_for_sorting.argsort(kind='mergesort')
        df = df.iloc[a]
        logger.info('Performed fancy sorting of chromosomes.')
    
    logger.info('Read %d lines (in %d chunks).', num_lines, num_chunks)
    logger.info('Found %d valid protein-coding gene entries.', c)
    logger.info('Final number of unique protein-coding genes: %d', df.shape[0])
    logger.info('Parsing time: %.1f s', t1-t0)
    
    # additional statistics
    all_chromosomes = list(df['Chromosome'].unique())
    logger.info('Valid chromosomes (%d): %s',
                len(all_chromosomes),
                ', '.join(all_chromosomes))
    logger.info('Excluded chromosomes (%d): %s',
                len(excluded_chromosomes),
                ', '.join(sorted(excluded_chromosomes)))
    
    logger.info('Sources:')
    for i, c in df['Source'].value_counts().iteritems():
        logger.info('\t%s: %d', i, c)
        
    logger.info('Gene types:')
    for i, c in df['Type'].value_counts().iteritems():
        logger.info('\t%s: %d', i, c)
    
    return df
