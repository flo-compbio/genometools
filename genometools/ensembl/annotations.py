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

# import os
import ftplib
import time
import re
import logging
from collections import Iterable, OrderedDict

import pandas as pd
    
# from .. import misc
from .. import gtf
from . import util

logger = logging.getLogger(__name__)


def get_annotation_urls_and_checksums(species, release=None, ftp=None):
    """Get FTP URLs and checksums for Ensembl genome annotations.
    
    Parameters
    ----------
    species : str or list of str
        The species or list of species for which to get genome annotations
        (e.g., "Homo_sapiens").
    release : int, optional
        The release number to look up. If `None`, use latest release. [None]
    ftp : ftplib.FTP, optional
        The FTP connection to use. If `None`, the function will open and close
        its own connection using user "anonymous".
    """
    ### type checks
    assert isinstance(species, (str, _oldstr)) or isinstance(species, Iterable)
    if release is not None:
        assert isinstance(release, int)
    if ftp is not None:
        assert isinstance(ftp, ftplib.FTP)

    ### open FTP connection if necessary
    close_connection = False
    ftp_server = 'ftp.ensembl.org'
    ftp_user = 'anonymous'
    if ftp is None:
        ftp = ftplib.FTP(ftp_server)
        ftp.login(ftp_user)
        close_connection = True    

    ### determine release if necessary
    if release is None:
        # use latest release
        release = util.get_latest_release(ftp=ftp)

    species_data = OrderedDict()
    if isinstance(species, (str, _oldstr)):
        species_list = [species]
    else:
        species_list = species
    for spec in species_list:

        # get the GTF file URL
        # => since the naming scheme isn't consistent across species,
        #    we're using a flexible scheme here to find the right file
        species_dir = '/pub/release-%d/gtf/%s' % (release, spec.lower())
        data = []
        ftp.dir(species_dir, data.append)
        gtf_file = []
        for d in data:
            i = d.rindex(' ')
            fn = d[(i + 1):]
            if fn.endswith('.%d.gtf.gz' % release):
                gtf_file.append(fn)
        assert len(gtf_file) == 1
        gtf_file = gtf_file[0]
        logger.debug('GTF file: %s', gtf_file)

        ### get the checksum for the GTF file
        checksum_url = '/'.join([species_dir, 'CHECKSUMS'])
        file_checksums = util.get_file_checksums(checksum_url, ftp=ftp)
        gtf_checksum = file_checksums[gtf_file]
        logger.debug('GTF file checksum: %d', gtf_checksum)

        gtf_url = 'ftp://%s%s/%s' %(ftp_server, species_dir, gtf_file)

        species_data[spec] = (gtf_url, gtf_checksum)

    # close FTP connection, if we opened it
    if close_connection:
        ftp.close()

    return species_data


def get_protein_coding_genes(
        path_or_buffer, chunksize=10000,
        chromosome_pattern=None,
        #chromosome_pattern=r'(?:\d\d?|MT|X|Y)$',
        include_polymorphic_pseudogenes=True,
        only_manual=False,
        remove_duplicates=True,
        sort_by='name'):
    r"""Get list of all protein-coding genes based on Ensembl GTF file.
    
    Parameters
    ----------
    path_or_buffer : str or buffer
        The GTF file (either the file path or a buffer)
    chromosome_pattern : str, optional
        Regular expression specifying valid chromosomes. [None]
    include_polymorphic_pseudogene : bool, optional
        Whether to include genes annotated as "polymorphic pseudogenes"?
    only_manual : bool, optional
        Whether to exclude annotations with source "ensembl", which
        are based only on an automatic annotation pipeline. [True]
    remove_duplicates : bool, optional
        Whether to remove duplicate annotations, i.e. those with different
        Ensembl IDs for the same gene. [True]
    sort_by : str, optional
        How to sort the genes. One of:
          - 'name': Genes are ordered alphabetically by their name
          - 'position': Genes are sorted by their position in the genome.abs
                        Genes are first sorted by chromosome, then by their
                        starting base pair position on the chromosome.
          - 'position_fancy': Like 'positional', but attempts to sort the
                              chromosomes in a more logical order than strictly
                              alphabetically. This currently works for human
                              and mouse genomes.abs
          - 'none': The order from the GTF file is retained. 
        Default: 'name'  

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
    chrompat = None
    if chromosome_pattern is not None:
        chrompat = re.compile(chromosome_pattern)
    
    c = 0
    num_lines = 0
    num_chunks = 0
    
    t0 = time.time()
    reader = pd.read_csv(path_or_buffer, encoding='ascii', sep='\t',
                         header=None, comment='#', dtype={0: str},
                         chunksize=chunksize)
    data = []
    header = ['name', 'ensembl_id',
              'chromosome', 'position', 'length',
              'source', 'type']
    
    valid_biotypes = set(['protein_coding'])
    if include_polymorphic_pseudogenes:
        valid_biotypes.add('polymorphic_pseudogene')
        
    valid_sources = set(['ensembl_havana', 'havana', 'insdc'])
    if not only_manual:
        # we also accept annotations with source "ensembl", which are the
        # product of an autmated annotation pipeline
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
            if chrompat is not None:
                match = chrompat.match(chrom)
                if match is None:
                    excluded_chromosomes.add(chrom)
                    continue

            c += 1

            ensembl_id = attr['gene_id']
            try:
                gene_name = attr['gene_name']
            except KeyError as e:
                gene_name = ensembl_id

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
        # make sure we only keep annotations with source "ensembl"
        # if no manual annotations are available
        sel = df['source'] == 'ensembl'
        redundant_ensembl_genes = set(df.loc[sel, 'name'].values) & \
                set(df.loc[~sel, 'name'].values)
        sel = sel & df['name'].isin(redundant_ensembl_genes)
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

        # sort by signed position value,
        # in order to make sure we keep the most "upstream" annotation in
        # the next step
        df.sort_values('position', kind='mergesort', inplace=True)

        # remove duplicates by keeping the first occurrence
        #df.drop_duplicates(['chromosome', 'name'], inplace=True)
        df.drop_duplicates('name', inplace=True)

        # restore original order using the numeric index
        df.sort_index(inplace=True)

        num_genes_after = df.shape[0]
        logger.info('Removed %d duplicate gene entries',
                    num_genes_before-num_genes_after)

    if sort_by == 'name':
        # sort alphabetically by gene name
        df.sort_values(['name'], kind='mergesort', inplace=True)

    elif sort_by in ['position', 'position_fancy']:
        # sort first by chromsome, then by absolute position
        df_sort = pd.concat([df['chromosome'], df['position'].abs()], axis=1)
        df_sort = df_sort.sort_values(['chromosome', 'position'],
                                      kind='mergesort')
        df = df.loc[df_sort.index]

        if sort_by == 'position_fancy':
            # Perform "fancy" positional sorting. Numbered chromosomes
            # are ordered numerically, and followed by the X, Y, and MT
            # chromosomes.
            def transform_chrom(chrom):
                """Helper function to obtain specific sort order."""
                try:
                    c = int(chrom)
                except:
                    if chrom in ['X', 'Y']:
                        return chrom
                    elif chrom == 'MT':
                        return '_MT'  # sort to the end
                    else:
                        return '__' + chrom  # sort to the very end
                else:
                    # make sure numbered chromosomes are sorted numerically
                    return '%02d' % c

            chrom_for_sorting = df['chromosome'].apply(transform_chrom)
            a = chrom_for_sorting.argsort(kind='mergesort')
            df = df.iloc[a]
            logger.info('Performed fancy sorting of chromosomes.')
    
    logger.info('Read %d lines (in %d chunks).', num_lines, num_chunks)
    logger.info('Found %d valid protein-coding gene entries.', c)
    logger.info('Final number of unique protein-coding genes: %d', df.shape[0])
    logger.info('Parsing time: %.1f s', t1-t0)
    
    # additional statistics
    all_chromosomes = list(df['chromosome'].unique())
    logger.info('Valid chromosomes (%d): %s',
                len(all_chromosomes),
                ', '.join(all_chromosomes))
    logger.info('Excluded chromosomes (%d): %s',
                len(excluded_chromosomes),
                ', '.join(sorted(excluded_chromosomes)))
    
    logger.info('Sources:')
    for i, c in df['source'].value_counts().iteritems():
        logger.info('\t%s: %d', i, c)
        
    logger.info('Gene types:')
    for i, c in df['type'].value_counts().iteritems():
        logger.info('\t%s: %d', i, c)
    
    return df
