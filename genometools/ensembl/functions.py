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
import sys
import ftplib
import re
import logging
from collections import Iterable, OrderedDict

from .. import misc

logger = logging.getLogger(__file__)

#ftp_server = 'ftp.ensembl.org'
#user = 'anonymous'


def get_latest_release(ftp=None):
    """Use files on the Ensembl FTP server to determine the latest release.

    Parameters
    ----------
    ftp : ftplib.FTP, optional
        FTP connection (with logged in user "anonymous").

    Returns
    -------
    int
        The version number of the latest release.
    """
    if ftp is not None:
        assert isinstance(ftp, ftplib.FTP)

    close_connection = False
    if ftp is None:
        ftp_server = 'ftp.ensembl.org'
        user = 'anonymous'
        ftp = ftplib.FTP(ftp_server)
        ftp.login(user)
        close_connection = True

    data = []
    ftp.dir('pub', data.append)
    pat = re.compile(r'.* current_README -> release-(\d+)/README$')
    latest = []
    for d in data:
        m = pat.match(d)
        if m is not None:
            latest.append(int(m.group(1)))

    assert len(latest) == 1, len(latest)
    latest = latest[0]

    if close_connection:
        ftp.close()

    return latest


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
        release = get_latest_release(ftp=ftp)

    species_data = OrderedDict()
    if isinstance(species, (str, _oldstr)):
        species_list = [species]
    else:
        species_list = species
    for spec in species_list:

        ### get the GTF file URL
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
        cs_path = '/'.join([species_dir, 'CHECKSUMS'])
        data = []
        ftp.retrbinary('RETR %s' % cs_path, data.append)
        # print(data)
        data = ''.join(d.decode('utf-8') for d in data).split('\n')[:-1]
        file_checksums = {}
        for d in data:
            file_name = d[(d.rindex(' ') + 1):]
            sum_ = int(d[:d.index(' ')])
            file_checksums[file_name] = sum_
        gtf_checksum = file_checksums[gtf_file]
        logger.debug('GTF file checksum: %d', gtf_checksum)

        gtf_url = 'ftp://%s%s/%s' %(ftp_server, species_dir, gtf_file)

        species_data[spec] = (gtf_url, gtf_checksum)

    if close_connection:
        ftp.close()

    return species_data


def get_cdna_url(species, release=None, ftp=None):
    """Returns the URL for a cDNA file hosted on the Ensembl FTP server.

    Parameters
    ----------
    species: str
        The scientific name of the species. It should be all lower-case,
        and the genus and species parts should be separated by an underscore
        (e.g., "homo_sapiens").
    release: int or ``None``, optional
        The Ensembl release number. If ``None``, the latest release is used.
        [None]
    ftp: ftplib.FTP or ``None``, optional
        The FTP connection. If ``None``, create a new connection. [None]
    """    
    #species_list, release=None, ftp=None

    ### type checks
    assert isinstance(species, (str, _oldstr))
    if release is not None:
        assert isinstance(release, int)
    if ftp is not None:
        assert isinstance(ftp, ftplib.FTP)

    ### open FTP connection if necessary
    close_connection = False
    if ftp is None:
        ftp_server = 'ftp.ensembl.org'
        ftp_user = 'anonymous'
        ftp = ftplib.FTP(ftp_server)
        ftp.login(ftp_user)
        close_connection = True

    ### determine release if necessary
    if release is None:
        # use latest release
        release = get_latest_release(ftp=ftp)
    
    # check if species exists
    #species_dir = 
    fasta_dir = '/pub/release-%d/fasta' % release 
    ftp.cwd(fasta_dir)
    if not species in ftp.nlst():
        logger.error('Species "%s" not found on Ensembl FTP server.', species)
        fasta_url = 'ftp://%s%s' %(ftp_server, fasta_dir)
        raise ValueError('Species "%s" not found. '
                         'See %s for a list of species available.'
                         % (species, fasta_url))
    
    cdna_dir = '/pub/release-%d/fasta/%s/cdna' %(release, species)
    ftp.cwd(cdna_dir)
    files = ftp.nlst()
    cdna_file = [f for f in files if f.endswith('.cdna.all.fa.gz')][0]

    cdna_url = 'ftp://%s%s/%s' %(ftp_server, cdna_dir, cdna_file)

    return cdna_url