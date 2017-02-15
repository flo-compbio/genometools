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

"""Utility functions for the Ensembl package."""

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
_oldstr = str
from builtins import *

import logging
import textwrap
import ftplib
import re
from collections import OrderedDict

from genometools import cli

logger = logging.getLogger(__name__)


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


def get_file_checksums(url, ftp=None):
    """Download and parse an Ensembl CHECKSUMS file and obtain checksums.

    Parameters
    ----------
    url : str
        The URL of the CHECKSUM file.
    ftp : `ftplib.FTP` or `None`, optional
        An FTP connection.
    
    Returns
    -------
    `collections.OrderedDict`
        An ordered dictionary containing file names as keys and checksums as
        values.

    Notes
    -----
    The checksums contains in Ensembl CHECKSUM files are obtained with the
    UNIX `sum` command.
    """
    assert isinstance(url, (str, _oldstr))
    if ftp is not None:
        assert isinstance(ftp, ftplib.FTP)

    # open FTP connection if necessary
    close_connection = False
    ftp_server = 'ftp.ensembl.org'
    ftp_user = 'anonymous'
    if ftp is None:
        ftp = ftplib.FTP(ftp_server)
        ftp.login(ftp_user)
        close_connection = True    
    
    # download and parse CHECKSUM file
    data = []
    ftp.retrbinary('RETR %s' % url, data.append)
    data = ''.join(d.decode('utf-8') for d in data).split('\n')[:-1]
    file_checksums = OrderedDict()
    for d in data:
        file_name = d[(d.rindex(' ') + 1):]
        sum_ = int(d[:d.index(' ')])
        file_checksums[file_name] = sum_
    
    logger.debug('Obtained checksums for %d files', len(file_checksums))

    # close FTP connection if we opened it
    if close_connection:
        ftp.close()
    
    return file_checksums


