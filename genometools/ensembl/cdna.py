
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

"""Functions for working with Ensembl cDNA data."""

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
_oldstr = str
from builtins import *

import logging
import ftplib

from . import util

logger = logging.getLogger(__name__)


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

    # type checks
    assert isinstance(species, (str, _oldstr))
    if release is not None:
        assert isinstance(release, int)
    if ftp is not None:
        assert isinstance(ftp, ftplib.FTP)

    # open FTP connection, if necessary
    close_connection = False
    if ftp is None:
        ftp_server = 'ftp.ensembl.org'
        ftp_user = 'anonymous'
        ftp = ftplib.FTP(ftp_server)
        ftp.login(ftp_user)
        close_connection = True

    # determine latest release, if necessary
    if release is None:
        release = util.get_latest_release(ftp=ftp)
    
    # check if species exists
    fasta_dir = '/pub/release-%d/fasta' % release 
    ftp.cwd(fasta_dir)
    if not species in ftp.nlst():
        logger.error('Species "%s" not found on Ensembl FTP server.', species)
        fasta_url = 'ftp://%s%s' %(ftp_server, fasta_dir)
        raise ValueError('Species "%s" not found. '
                         'See %s for a list of species available.'
                         % (species, fasta_url))

    # determine URL of the cdna file
    # (file names are not consistent across species)
    cdna_dir = '/pub/release-%d/fasta/%s/cdna' %(release, species)
    ftp.cwd(cdna_dir)
    files = ftp.nlst()
    cdna_file = [f for f in files if f.endswith('.cdna.all.fa.gz')][0]
    cdna_url = 'ftp://%s%s/%s' %(ftp_server, cdna_dir, cdna_file)

    # close FTP connection, if we opened it
    if close_connection:
        ftp.close()

    return cdna_url
