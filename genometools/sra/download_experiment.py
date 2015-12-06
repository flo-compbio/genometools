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

"""Module containing the ``sra_download_experiment.py`` script.

"""

import sys
import os
import argparse
import logging
import ftputil

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
    parser = argparse.ArgumentParser(description=
        'Downloads all .sra files from NCBI SRA for a given experiment ID.')

    parser.add_argument('-e','--experiment-id', required=True,
        help='The SRA experiment ID (starting with "SRX").')

    parser.add_argument('-d','--download-dir', required=True,
        help="""Path of download directory. (The script will create a
                subdirectory in that directory named after the experiment
                ID.""")

    parser.add_argument('-o','--overwrite', action='store_true',
        help="""Allow files to be overwritten. If not specified, a warning will
                be issued for each existing file, and the file will not be
                downloaded.""")

    parser.add_argument('-l','--log-file', default=None,
        help='Path of log file. If not specified, print to stdout.')

    parser.add_argument('-q','--quiet', action='store_true',
        help='Suppress all output except warnings and errors.')

    parser.add_argument('-v','--verbose', action='store_true',
        help='Enable verbose output. Ignored if ``--quiet`` is specified.')

    return parser

def sizeof_fmt(num, suffix='B'):
    for unit in ['','Ki','Mi','Gi','Ti','Pi','Ei','Zi']:
        if abs(num) < 1024.0:
            return "%3.0f%s%s" % (num, unit, suffix)
        num /= 1024.0
    return "%.0f%s%s" % (num, 'Yi', suffix)

def main(args=None):
    """Download all .sra from NCBI SRA for a given experiment ID.

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
    if args is None:
        # parse command-line arguments
        parser = get_argument_parser()
        args = parser.parse_args()

    experiment_id = args.experiment_id
    download_dir = args.download_dir
    overwrite = args.overwrite
    log_file = args.log_file
    quiet = args.quiet
    verbose = args.verbose

    # configure logger
    log_stream = sys.stdout

    log_level = logging.INFO
    if quiet:
        log_level = logging.WARNING
    elif verbose:
        log_level = logging.DEBUG
    logger = misc.configure_logger(__name__, log_stream = log_stream,
            log_file = log_file, log_level = log_level)

    host = 'ftp-trace.ncbi.nlm.nih.gov'
    user = 'anonymous'
    password = 'anonymous'

    output_dir = download_dir + experiment_id + '/'
    # make sure output directory exists
    misc.make_sure_dir_exists(output_dir)
    logger.info('Created output directory: "%s".', output_dir)

    with ftputil.FTPHost(host,user,password) as ftp_host:
        download_paths = []
        exp_dir = '/sra/sra-instant/reads/ByExp/sra/SRX/%s/%s/' \
                %(experiment_id[:6],experiment_id)
        ftp_host.chdir(exp_dir)
        run_folders = ftp_host.listdir(ftp_host.curdir)
        logging.info('Found %d run folders.',len(run_folders))

        # compile a list of all files to be downloaded
        for folder in run_folders:
            files = ftp_host.listdir(folder)
            assert len(files) == 1
            for f in files:
                download_paths.append([folder,f])

        # download files
        n_files = len(download_paths)
        logger.info('Found %d files to be downloaded.', n_files)
        for i,path in enumerate(download_paths):
            logger.debug('%02d - %s/%s' %(i+1,path[0],path[1]))

        c = 0
        #def keepalive(chunk):
        #    c += chunk
        #    print '\r%d' %(c
        #    ftp_host.keep_alive()

        for i,path in enumerate(download_paths):
            file_name = path[1]
            source = path[0] + '/' + file_name
            output_file = output_dir + file_name
            # get file size
            stat = ftp_host.stat(source)
            size_str = sizeof_fmt(stat.st_size)
            c = 0
            if os.path.isfile(output_file) and not overwrite:
                logger.warning('Output file "%s" already exists! Skipping...',
                        output_file)
            else:
                logger.info('Downloading file %d / %d: %s (%s)', i+1, n_files,
                        file_name, size_str)
                ftp_host.download(source,output_file)
                #with ftp_host.open(source, 'rb') as source, \
                #        open(output_file,'wb') as target:
                #    ftp_host.copyfileobj(source,target)
                

    return 0

if __name__ == '__main__':
    return_code = main()
    sys.exit(return_code)
