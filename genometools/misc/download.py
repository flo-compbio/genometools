# Copyright (c) 2015-2017 Florian Wagner
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

"""HTTP and FTP download functions."""

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
_oldstr = str
from builtins import *

import logging
import os
import shutil
import ftplib

import six
if six.PY3:
    from urllib import parse as urlparse
else:
    import urlparse

import requests

_logger = logging.getLogger(__name__)


def http_download(url, download_file,
                  overwrite=False, raise_http_exception=True):
    """Download a file over HTTP(S).
    
    See: http://stackoverflow.com/a/13137873/5651021 

    Parameters
    ----------
    url : str
        The URL.
    download_file : str
        The path of the local file to write to.
    overwrite : bool, optional
        Whether to overwrite an existing file (if present). [False]
    raise_http_exception : bool, optional
        Whether to raise an exception if there is an HTTP error. [True]

    Raises
    ------
    OSError
        If the file already exists and overwrite is set to False.
    `requests.HTTPError`
        If an HTTP error occurred and `raise_http_exception` was set to `True`.
    """

    assert isinstance(url, (str, _oldstr))
    assert isinstance(download_file, (str, _oldstr))
    assert isinstance(overwrite, bool)
    assert isinstance(raise_http_exception, bool)

    u = urlparse.urlparse(url)
    assert u.scheme in ['http', 'https']

    if os.path.isfile(download_file) and not overwrite:
        raise OSError('File "%s" already exists!' % download_file)
    
    r = requests.get(url, stream=True)
    if raise_http_exception:
        r.raise_for_status()
    if r.status_code == 200:
        with open(download_file, 'wb') as fh:
            r.raw.decode_content = True
            shutil.copyfileobj(r.raw, fh)
        _logger.info('Downloaded file "%s".', download_file)
    else:
        _logger.error('Failed to download url "%s": HTTP status %d/%s',
                      url, r.status_code, r.reason)


def ftp_download(url, download_file, if_exists='error',
                 user_name='anonymous', password='', blocksize=4194304):
    """Downloads a file from an FTP server.

    Parameters
    ----------
    url : str
        The URL of the file to download.
    download_file : str
        The path of the local file to download to. 
    if_exists : str, optional
        Desired behavior when the download file already exists. One of:
          'error'     - Raise an OSError
          'skip'      - Do nothing, only report a warning.
          'overwrite' - Overwrite the file. reporting a warning.
        Default: 'error'.
    user_name : str, optional
        The user name to use for logging into the FTP server. ['anonymous']
    password : str, optional
        The password to use for logging into the FTP server. ['']
    blocksize : int, optional
        The blocksize (in bytes) to use for downloading. [4194304]

    Returns
    -------
    None
    """
    assert isinstance(url, (str, _oldstr))
    assert isinstance(download_file, (str, _oldstr))
    assert isinstance(if_exists, (str, _oldstr))
    assert isinstance(user_name, (str, _oldstr))
    assert isinstance(password, (str, _oldstr))

    u = urlparse.urlparse(url)
    assert u.scheme == 'ftp'

    if if_exists not in ['error', 'skip', 'overwrite']:
        raise ValueError('"if_exists" must be "error", "skip", or "overwrite" '
                         '(was: "%s").', str(if_exists))

    if os.path.isfile(download_file):
        if if_exists == 'error':
            raise OSError('File "%s" already exists.' % download_file)
        elif if_exists == 'skip':
            _logger.warning('File "%s" already exists. Skipping...',
                            download_file)
            return
        else:
            _logger.warning('Overwriting file "%s"...', download_file)

    ftp_server = u.netloc
    ftp_path = u.path

    if six.PY3:
        with ftplib.FTP(ftp_server) as ftp:
            ftp.login(user_name, password)
            with open(download_file, 'wb') as ofh:
                ftp.retrbinary('RETR %s' % ftp_path,
                               callback=ofh.write, blocksize=blocksize)
    else:
        ftp = ftplib.FTP(ftp_server)
        ftp.login(user_name, password)
        with open(download_file, 'wb') as ofh:
            ftp.retrbinary('RETR %s' % ftp_path,
                           callback=ofh.write, blocksize=blocksize)
        ftp.close()

    _logger.info('Downloaded file "%s" over FTP.', download_file)
