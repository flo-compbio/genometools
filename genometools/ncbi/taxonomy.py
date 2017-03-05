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

"""Functions for parsing data from the NCBI Taxonomy database."""

import tarfile
from collections import defaultdict

import pandas as pd


def _get_divisions(taxdump_file):
    """Returns a dictionary mapping division names to division IDs."""
    
    with tarfile.open(taxdump_file) as tf:    
        with tf.extractfile('division.dmp') as fh:
            df = pd.read_csv(fh, header=None, sep='|', encoding='ascii')

    # only keep division ids and names
    df = df.iloc[:, [0, 2]]

    # remove tab characters flanking each division name
    df.iloc[:, 1] = df.iloc[:, 1].str.strip('\t')

    # generate dictionary
    divisions = {}
    for _, row in df.iterrows():
        divisions[row.iloc[1]] = row.iloc[0]
    
    return divisions


def _get_species_taxon_ids(taxdump_file,
                           select_divisions=None, exclude_divisions=None):
    """Get a list of species taxon IDs (allow filtering by division)."""
    
    if select_divisions and exclude_divisions:
        raise ValueError('Cannot specify "select_divisions" and '
                         '"exclude_divisions" at the same time.')

    select_division_ids = None
    exclude_division_ids = None
    
    divisions = None
    if select_divisions or exclude_divisions:
        divisions = _get_divisions(taxdump_file)
        
    if select_divisions:
        select_division_ids = set([divisions[d] for d in select_divisions])
        
    elif exclude_divisions:
        exclude_division_ids = set([divisions[d] for d in exclude_divisions])
    
    with tarfile.open(taxdump_file) as tf:    
    
        with tf.extractfile('nodes.dmp') as fh:
            df = pd.read_csv(fh, header=None, sep='|', encoding='ascii')

    # select only tax_id, rank, and division id columns
    df = df.iloc[:, [0, 2, 4]]

    if select_division_ids:
        # select only species from specified divisions
        df = df.loc[df.iloc[:, 2].isin(select_division_ids)]
            
    elif exclude_division_ids:
        # exclude species from specified divisions
        df = df.loc[~df.iloc[:, 2].isin(exclude_division_ids)]
    
    # remove tab characters flanking each rank name
    df.iloc[:, 1] = df.iloc[:, 1].str.strip('\t')

    # get taxon IDs for all species
    taxon_ids = df.iloc[:, 0].loc[df.iloc[:, 1] == 'species'].values
    return taxon_ids


#exclude_divisions = ['Bacteria', 'Phages', 'Synthetic and Chimeric',
#                     'Unassigned', 'Viruses', 'Environmental samples']

def get_species(taxdump_file, select_divisions=None,
                exclude_divisions=None, nrows=None):
    """Get a dataframe with species information."""
    
    if select_divisions and exclude_divisions:
        raise ValueError('Cannot specify "select_divisions" and '
                         '"exclude_divisions" at the same time.')

    select_taxon_ids = _get_species_taxon_ids(
        taxdump_file,
        select_divisions=select_divisions,
        exclude_divisions=exclude_divisions)
    select_taxon_ids = set(select_taxon_ids)
    
    with tarfile.open(taxdump_file) as tf:
        with tf.extractfile('names.dmp') as fh:    
            df = pd.read_csv(fh, header=None, sep='|',
                             encoding='ascii', nrows=nrows)

    # only keep information we need
    df = df.iloc[:, [0, 1, 3]]

    # only select selected species
    df = df.loc[df.iloc[:, 0].isin(select_taxon_ids)]

    # remove tab characters flanking each "name class" entry
    df.iloc[:, 2] = df.iloc[:, 2].str.strip('\t')

    # select only "scientific name" and "common name" rows
    df = df.loc[df.iloc[:, 2].isin(['scientific name', 'common name'])]

    # remove tab characters flanking each "name" entry 
    df.iloc[:, 1] = df.iloc[:, 1].str.strip('\t')
    
    # collapse common names for each scientific name
    common_names = defaultdict(list)
    cn = df.loc[df.iloc[:, 2] == 'common name']
    for _, row in cn.iterrows():
        common_names[row.iloc[0]].append(row.iloc[1])
        
    # build final dataframe (this is very slow)
    sn = df.loc[df.iloc[:, 2] == 'scientific name']
    species = []
    for i, row in sn.iterrows():
        species.append([row.iloc[0], row.iloc[1],
                        '|'.join(common_names[row.iloc[0]])])
    species_df = pd.DataFrame(species).set_index(0)
    species_df.columns = ['scientific_name', 'common_names']
    species_df.index.name = 'taxon_id'
    return species_df
