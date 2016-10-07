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
from builtins import *

import os
import sys
import ftplib
import re
from collections import Iterable

from .. import misc

logger = misc.get_logger()

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


def download_gene_annotations(species, download_dir, release=None,
                              ftp=None, avoid_redownload=True):
    """Download gene annotations (in GTF format) from Ensembl FTP server.

    Parameters
    ----------
    species : list of str
        List of species for which to download gene annotations.
        Each species must be specified by its scientific name, with an
        undescore separating genus and species. (e.g., 'Homo_sapiens').
    download_dir : str
        The path of the download directory. Must be an existing directory.
    release : int, optional
        The Ensembl release number.
        If ``None``, download data for the latest release. [None]
    ftp : ftplib.FTP, optional
        The FTP connection to use. [None]
    avoid_redownload : bool
        Whether to avoid re-downloading files if their checksums agree.
        Only available on platforms that have the "sum" unix program.

    Returns
    -------
    list of str
        The paths of the downloaded file(s).
    """
    assert isinstance(species, Iterable)
    if release is not None:
        assert isinstance(release, int)
    assert isinstance(download_dir, str)

    if not os.path.isdir(download_dir):
        raise OSError('Download directory "%s" does not exist.' % download_dir)

    if ftp is not None:
        assert isinstance(ftp, ftplib.FTP)

    close_connection = False
    if ftp is None:
        ftp_server = 'ftp.ensembl.org'
        user = 'anonymous'
        ftp = ftplib.FTP(ftp_server)
        ftp.login(user)
        close_connection = True

    if release is None:
        # use latest release
        release = get_latest_release(ftp=ftp)

    # determine whether we can call the unix "sum" program to calculate checksums
    test_checksum = False
    if sys.platform.startswith('linux') or sys.platform in ['cygwin',
                                                            'darwin']:
        test_checksum = True

    downloaded_files = []

    ### update gene annotations for all species
    for spec in species:
        logger.info('Now processing species "%s"', spec)

        # spec_name = species_names[spec]

        # find the precise name of the GTF file that we're interested in
        # Note: Ensembl uses lower-case scientific species names as folder names
        #      e.g., "homo_sapiens").
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
        # self.logger.debug(gtf_file)

        # download the CHECKSUMS file and create a mapping of file names to checksums
        cs_path = '/'.join([species_dir, 'CHECKSUMS'])
        data = []
        ftp.retrbinary('RETR %s' % cs_path, data.append)
        print(data)
        data = ''.join(d.decode('utf-8') for d in data).split('\n')[:-1]
        checksums = {}
        for d in data:
            file_name = d[(d.rindex(' ') + 1):]
            sum_ = int(d[:d.index(' ')])
            checksums[file_name] = sum_
        print(checksums)

        # compare checksums to see if we need to download the file
        download_file = os.path.join(download_dir, gtf_file)
        if avoid_redownload and test_checksum and os.path.isfile(download_file) \
                and misc.test_file_checksum(download_file,
                                            checksums[gtf_file]):
            logger.info(
                'File "%s" already present and checksums agree. Skipping download.',
                gtf_file)
        else:
            download_path = '/'.join([species_dir, gtf_file])
            logger.debug('Download path: %s', download_path)
            logger.debug('Downloading to file: %s', download_file)
            with open(download_file, 'wb') as ofh:
                ftp.retrbinary('RETR %s' % download_path, ofh.write)

            logger.debug('Done!')
            success = True
            if test_checksum:
                if not misc.test_file_checksum(download_file,
                                               checksums[gtf_file]):
                    logger.error(
                        'Checksums don''t agree! Deleting downloaded file...')
                    success = False
            else:
                if (not os.path.isfile(download_file)) or os.path.getsize(
                        download_file) == 0:
                    success = False

            if not success:
                logger.error('Something went wrong with the file download.')
                try:
                    os.remove(download_file)
                except OSError:
                    pass
                download_file = None
        if download_file is not None:
            downloaded_files.append(download_file)

    if close_connection:
        ftp.close()

    return downloaded_files