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

import sys
import os
import logging
import argparse
import csv
import gzip

def get_argument_parser():
    parser = argparse.ArgumentParser(description='Generate a mapping of Entrez IDs to gene symbols.')

    parser.add_argument('-f','--gene2acc-file',required=True,'gene2accession.gz file from ftp://ftp.ncbi.nlm.nih.gov/gene/DATA, or a subset thereof')
    parser.add_argument('-o','--output-file',required=True,'Output file')
    parser.add_argument('-l','--log-file',default=None,help='Log file - if not specified, print to stdout')
    parser.add_argument('-v','--verbose',action='store_true',help='Verbose output')

    return parser

def read_args_from_cmdline():
    parser = get_argument_parser()
    return parser.parse_args()

def open_plain_or_gzip(fn):
    try:
        gzip.open(fn).next()
        return gzip.open(fn)
    except IOError:
        return open(fn)

def read_gene2acc(fn,logger):
    gene2acc = {}
    with open_plain_or_gzip(fn) as fh:
        reader = csv.reader(fh,dialect='excel-tab')
        reader.next() # skip header
        for i,l in enumerate(reader):
            if (i % 1000000) == 0:
                logger.info('%d...',i)

            id_ = int(l[1])
            symbol = l[15]

            try:
                gene2acc[id_].append(symbol)
            except:
                gene2acc[id_] = [symbol]
            #print l[0],l[15]
    logger.info('')

    # make sure all EntrezIDs map to a unique gene symbol
    n = len(gene2acc.keys())
    for k,v in gene2acc.iteritems():
        symbols = sorted(set(v))
        assert len(symbols) == 1
        gene2acc[k] = symbols[0]

    all_symbols = sorted(set(gene2acc.values()))
    m = len(all_symbols)

    logger.info('Found %d Entrez Gene IDs associated with %d gene symbols.', n,m)
    return gene2acc

def write_entrez2gene(ofn,entrez2gene,logger):
    with open(ofn,'w') as ofh:
        writer = csv.writer(ofh,dialect='excel-tab',lineterminator='\n')
        for k,v in entrez2gene.iteritems():
            writer.writerow([k,v])

def main(args=None):

    if args is None:
        args = read_args_from_cmdline()

    gene2acc_file = args.gene2acc_file
    output_file = args.output_file
    log_file = args.log_file
    verbose = args.verbose

    # configure logger
    log_level = logging.INFO
    if verbose:
        log_level = logging.DEBUG

    log_format = '[%(asctime)s] %(levelname)s: %(message)s'
    log_datefmt = '%Y-%m-%d %H:%M:%S'
    # when filename is not None, "stream" parameter is ignored (see https://docs.python.org/2/library/logging.html#logging.basicConfig)
    logging.basicConfig(filename=log_file,stream=sys.stdout,level=log_level,format=log_format,datefmt=log_datefmt)
    logger = logging.getLogger()

    entrez2gene = read_gene2acc(gene2acc_file,logger)
    write_entrez2gene(output_file,entrez2gene,logger)

    return 0

if __name__ == '__main__':
    return_code = main()
    sys.exit(return_code)
