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

"""Functions for downloading files from GDC."""

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
_oldstr = str
from builtins import *

import os
import logging
import hashlib
from contextlib import closing
from collections import OrderedDict

import pandas as pd
import requests

logger = logging.getLogger(__name__)


def download_files(manifest, download_dir, auth_token=None,
                       chunk_size=1048576, avoid_redownload=True):
    """Individually download files from GDC.
    
    Params
    ------
    manifest : `pandas.DataFrame`
        GDC manifest that contains a list of files. The data frame should
        have five columns: id, filename, md5, size, and state.
    download_dir : str
        The path of the download directory.
    auth_token : str, optional
        Authentication token for downloading protected data.
        If None, do not send authentication header. [None]
    chunk_size : int, optional
        The chunk size (in bytes) to use for downloading data. [1048576]
        
    Returns
    -------
    None
    """
    assert isinstance(manifest, pd.DataFrame)
    assert isinstance(download_dir, str)
    assert isinstance(chunk_size, int)
    if auth_token is not None:
        assert isinstance(auth_token, str)
                       
    def get_file_md5hash(path):
        """Calculate the MD5 hash for a file."""
        with open(path, 'rb') as fh:
            h = hashlib.md5(fh.read()).hexdigest()
        return h

    headers = {}
    if auth_token is not None:
        headers['X-Auth-Token'] = auth_token
        
    #payload = {'ids': file_ids}    
    #logger.info('Downloading data to "%s"...', download_dir)
    num_files = manifest.shape[0]
    logger.info('Downloading %d files to "%s"...', num_files, download_dir)
    for i, row in manifest.iterrows():
        #(uuid, (file_name, file_hash)) 
        success = False        
        download_file = os.path.join(download_dir, row['filename'])
        if ((i+1) % 100) == 0:
            logger.info('Downloading file %d / %d...', i+1, num_files)
        
        if avoid_redownload and os.path.isfile(download_file) and \
                get_file_md5hash(download_file) == row['md5']:
            logger.info('File %s already downloaded...skipping.',
                        download_file)
            success = True
        
        while not success:
            with closing(
                requests.get('https://gdc-api.nci.nih.gov/data/%s'
                             % row['id'], headers=headers, stream=True)) \
                as r:
                # get suggested file name from "Content-Disposition" header
                # suggested_file_name = re.findall(
                #       "filename=(.+)", r.headers['Content-Disposition'])[0]

                r.raise_for_status()

                with open(download_file, 'wb') as ofh:
                    for chunk in r.iter_content(chunk_size=chunk_size): 
                        if chunk: # filter out keep-alive new chunks
                            ofh.write(chunk)
            with open(download_file, 'rb') as fh:
                h = hashlib.md5(fh.read()).hexdigest()
                if h == row['md5']:
                    success = True
            if not success:
                logger.warning('Hash value mismatch (should be: %s; is: %s). '
                               'Attempting to re-download file...',
                               row['md5'], h)
