#!/usr/bin/env python3
#cython: boundscheck=False

import numpy as np
cimport numpy as np
np.import_array()

import sys
import os
import logging
import argparse
import subprocess
import io
import itertools as it
import gzip
import errno

import pandas as pd
import numpy as np


def make_sure_dir_exists(dir_, create_subfolders=False):
    """Ensures that a directory exists.

    Adapted from StackOverflow users "Bengt" and "Heikki Toivonen"
    (http://stackoverflow.com/a/5032238).

    Parameters
    ----------
    dir_: str
        The directory path.
    create_subfolders: bool, optional
        Whether to create any inexistent subfolders. [False]
    
    Returns
    -------
    None

    Raises
    ------
    OSError
        If a file system error occurs.
    """
    assert isinstance(dir_, str)
    assert isinstance(create_subfolders, bool)

    try:
        if create_subfolders:
            os.makedirs(dir_)
        else:
            os.mkdir(dir_)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise


def configure_logger(name, log_stream=sys.stdout, log_file=None,
        log_level=logging.INFO, keep_old_handlers=False,
        propagate=False):
    """Configures and returns a logger.

    This function serves to simplify the configuration of a logger that
    writes to a file and/or to a stream (e.g., stdout).

    Parameters
    ----------
    name: str
        The name of the logger. Typically set to ``__name__``.
    log_stream: a stream object, optional
        The stream to write log messages to. If ``None``, do not write to any
        stream. The default value is `sys.stdout`.
    log_file: str, optional
        The path of a file to write log messages to. If None, do not write to
        any file. The default value is ``None``.
    log_level: int, optional
        A logging level as `defined`__ in Python's logging module. The default
        value is `logging.INFO`.
    keep_old_handlers: bool, optional
        If set to ``True``, keep any pre-existing handlers that are attached to
        the logger. The default value is ``False``.
    propagate: bool, optional
        If set to ``True``, propagate the loggers messages to the parent logger.
        The default value is ``False``.

    Returns
    -------
    `logging.Logger`
        The logger.

    Notes
    -----
    Note that if ``log_stream`` and ``log_file`` are both ``None``, no handlers
    will be created.

    __ loglvl_

    .. _loglvl: https://docs.python.org/2/library/logging.html#logging-levels

    """
    # create a child logger
    logger = logging.getLogger(name)

    # set the logger's level
    logger.setLevel(log_level)

    # set the logger's propagation attribute
    logger.propagate = propagate

    if not keep_old_handlers:
        # remove previously attached handlers
        logger.handlers = []

    # create the formatter
    log_fmt = '[%(asctime)s] %(levelname)s: %(message)s'
    log_datefmt = '%Y-%m-%d %H:%M:%S'
    formatter = logging.Formatter(log_fmt,log_datefmt)

    # create and attach the handlers

    if log_stream is not None:
        # create a StreamHandler
        stream_handler = logging.StreamHandler(log_stream)
        stream_handler.setFormatter(formatter)
        logger.addHandler(stream_handler)

    if log_file is not None:
        # create a FileHandler
        file_handler = logging.FileHandler(log_file)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    if log_stream is None and log_file is None:
        # "no handler" => use NullHandler
        logger.addHandler(logging.NullHandler())
      
    return logger


def get_logger(name='', log_stream=None, log_file=None,
        quiet=False, verbose=False):
    """Convenience function for getting a logger."""

    # configure root logger
    log_level = logging.INFO
    if quiet:
        log_level = logging.WARNING
    elif verbose:
        log_level = logging.DEBUG

    if log_stream is None:
        log_stream = sys.stdout

    new_logger = configure_logger(name, log_stream=log_stream,
            log_file=log_file, log_level=log_level)

    return new_logger


def zcat_subproc(path):
    subproc = subprocess.Popen('gunzip -c "%s"' % path, shell=True,
                               stdout=subprocess.PIPE)
    return subproc


def get_reverse_complement(seq):
    rc = {
        'A': 'T',
        'T': 'A',
        'G': 'C',
        'C': 'G'
    }
    compseq = ''.join([rc[nuc] for nuc in seq[::-1]])
    return compseq
        

def get_mismatch_sequences(seq):
    # get all nucleotide sequences with hamming distance 1 to seq
    for pos in range(len(seq)):
        for nuc in ['A', 'C', 'G', 'T']:
            if nuc == seq[pos]:
                continue
            mm = seq[:pos] + nuc + seq[(pos+1):]
            yield mm


def get_barcode_exact_mapping(barcode_file):
    """
    Returns a mapping from sequences to barcodes, for
    all sequences matching barcodes exactly.
    
    Only the first 8 nucleotides of each barcode are taken into account.
    """ 
    import pandas as pd
    from collections import defaultdict, Counter

    barcodes = pd.read_csv(barcode_file, squeeze=True, header=None)

    # take reverse complement sequence
    barcodes = barcodes.apply(get_reverse_complement)
    
    mapping = dict([bc, i] for i, bc in enumerate(barcodes))
    #mapping = dict([bc, (bc, len(bc), i)] for i, bc in enumerate(barcodes))
    return mapping


def get_barcode_mismatch_mapping(barcode_file):
    """
    Returns a mapping from sequences to barcodes, for
    all sequences with one mismatch, excluding abiguous cases.
    
    Only the first 8 nucleotides of each barcode are taken into account.
    """
    
    import pandas as pd
    from collections import defaultdict, Counter

    barcodes = pd.read_csv(barcode_file, squeeze=True, header=None)

    # take reverse complement sequence
    barcodes = barcodes.apply(get_reverse_complement)
                       
    # exact matches should always be preferred to 1-mismatch matches
    exact = set(barcodes)

    assert len(exact) == len(barcodes)

    # generate a list of all 1-mismatch mathces
    # to do so, generate a dictionary of "mismatch sequence"
    #   => ["barcode sequence 1", ...]
    # this lets us map back to how many mismatches were generated for each
    # barcode
    
    mismatch = defaultdict(set)
    #print(len(barcodes))
    for i, bc in enumerate(barcodes):
        for pos in range(len(bc)):
            for nuc in ['A', 'C', 'G', 'T']:
                mm = bc[:pos] + nuc + bc[(pos+1):]
                if mm not in exact:
                    mismatch[mm].add((bc, len(bc), i))
                    #mismatch[mm].add(i)

    LOGGER.info('How many barcodes are assigned to each mismatch sequence?')
    LOGGER.info(str(Counter(len(v) for v in mismatch.values())))
    # print([(k,v) for k,v in mismatch.items() if len(v) == 8])
    
    bc_count = Counter()
    for v in mismatch.values():
        for bc in v:
            bc_count[bc[0]] += 1
    LOGGER.info('In how many mismatch sequences does each barcode occur?')
    LOGGER.info(str(Counter(bc_count.values())))
    
    # only keep mismatch mapping
    mismatch_unambiguous = dict(
        (mm, list(barcodes)[0][2])
        for mm, barcodes in mismatch.items() if len(barcodes) == 1)
    LOGGER.info('Found %d unambiguously mapping sequences.',
                len(mismatch_unambiguous))
    return mismatch_unambiguous


def get_all_kmers(k, kmer='', kmer_list=None):
    if kmer_list is None:
        kmer_list = []
    if len(kmer) == k:
        kmer_list.append(kmer)
    else:
        for nuc in ['A', 'C', 'G', 'T']:
            var = kmer + nuc
            get_all_kmers(k, var, kmer_list)
    if not kmer:
        return kmer_list


def count_barcodes_cython(barcode_read_file, mrna_read_file,
                          barcode1_file, barcode2_file,
                          output_dir, logger,
                          max_reads=None):
    
    #1000000
    #cdef int MAX_READS = max_reads
    if max_reads is None:
        max_reads = int(2e9)

    # create output directory, if necessary
    make_sure_dir_exists(output_dir)
    
    # define output files
    output_count_file = os.path.join(output_dir, 'barcode_counts.tsv')
    output_read_file = os.path.join(output_dir, 'processed_reads.fastq')
    log_file = os.path.join(output_dir, 'log.txt')

    # add file handler to the logger
    file_handler = logging.FileHandler(log_file)
    file_handler.setFormatter(logger.handlers[0].formatter)
    logger.addHandler(file_handler)

    w1 = 'GAGTGATTGCTTGTGACGCCTT'
    
    w1_mapping = {w1}
    w1_mapping |= set(get_mismatch_sequences(w1))
    
    bc1_exact_mapping = get_barcode_exact_mapping(barcode1_file)
    bc1_mismatch_mapping = get_barcode_mismatch_mapping(barcode1_file)
    
    bc2_exact_mapping = get_barcode_exact_mapping(barcode2_file)
    bc2_mismatch_mapping = get_barcode_mismatch_mapping(barcode2_file)   
    
    # first barcodes
    barcodes1 = pd.read_csv(
        barcode1_file, sep='\t', header=None, squeeze=True) \
            .apply(get_reverse_complement).values
    assert barcodes1.ndim == 1
    cdef unsigned long [:] bc1_counts = \
            np.zeros(barcodes1.size, dtype=np.uint64)
    
    # second barcodes
    barcodes2 = pd.read_csv(
        barcode2_file, sep='\t', header=None, squeeze=True) \
            .apply(get_reverse_complement).values
    assert barcodes2.ndim == 1
    cdef unsigned long [:] bc2_counts = \
            np.zeros(barcodes2.size, dtype=np.uint64)
    
    cdef unsigned long [:, ::1] bc_counts = \
            np.zeros((barcodes1.size, barcodes2.size), dtype=np.uint64)
        
    subproc1 = zcat_subproc(barcode_read_file)
    subproc2 = zcat_subproc(mrna_read_file)

    wrap1 = io.TextIOWrapper(subproc1.stdout, encoding='ascii')
    wrap2 = io.TextIOWrapper(subproc2.stdout, encoding='ascii')
    
    umi_numbers = dict([seq, ind] for ind, seq in enumerate(get_all_kmers(6)))    
    
    cdef int i, c, bc1, bc2, umi_num
    cdef int w1_len
    cdef int bc1_found, bc2_found
    cdef int passed_w1, passed_bc
    cdef int bc1_len, bc2_start, umi_start

    c = 0
    w1_len = len(w1)
    passed_w1 = 0
    passed_bc = 0
    
    with open(output_read_file, 'w', encoding='ascii') as out:
        for i in range(max_reads):
            # skip over title line in FASTQ

            try:
                next(wrap1)
                l1 = next(wrap2)
            except StopIteration:
                break

            r1 = next(wrap1)
            r2 = next(wrap2)

            c += 1
            if (c % 10000000) == 0:
                LOGGER.info('Processed %d reads...', c)

            bc1_len = 0
            if len(r1) >= 48:  # make sure read is long enough
                if r1[8:(8+w1_len)] in w1_mapping:
                    bc1_len = 8
                elif r1[9:(9+w1_len)] in w1_mapping:
                    bc1_len = 9
                elif r1[10:(10+w1_len)] in w1_mapping:
                    bc1_len = 10
                elif r1[11:(11+w1_len)] in w1_mapping:
                    bc1_len = 11

            if bc1_len == 0:
                # could not identify w1 sequence
                # => skip this read (i.e., the next 2 lines in each file)
                next(wrap1)
                next(wrap1)
                next(wrap2)
                next(wrap2)                
                continue

            passed_w1 += 1

            bc1_found = 1
            bc2_found = 1

            # try to identify first barcode
            try:
                bc1 = bc1_exact_mapping[r1[:bc1_len]]
            except KeyError:
                try:
                    bc1 = bc1_mismatch_mapping[r1[:bc1_len]]
                except KeyError:
                    bc1_found = 0

            bc2_start = bc1_len + w1_len

            # try to identify second barcode
            try:
                bc2 = bc2_exact_mapping[r1[bc2_start:(bc2_start+8)]]
            except KeyError:
                try:
                    bc2 = bc2_mismatch_mapping[r1[bc2_start:(bc2_start+8)]]
                except KeyError:
                    bc2_found = 0

            # extract UMI sequence
            umi_start = bc2_start + 8
            umi = r1[umi_start:(umi_start+6)]

            # read the next 2 lines in each file
            next(wrap1)
            next(wrap1)
            l3 = next(wrap2)
            l4 = next(wrap2)

            if bc1_found != 0 and bc2_found != 0 and 'N' not in umi:
                # the read is good, we identified both barcodes,
                # and no 'N' in UMI
                passed_bc += 1

                bc1_counts[bc1] += 1
                bc2_counts[bc2] += 1
                bc_counts[bc1, bc2] += 1

                umi_num = umi_numbers[umi]

                out.write('@%03d-%03d_%04d_%s' %(bc1, bc2, umi_num, l1[1:]))
                out.write(r2)
                out.write(l3)
                out.write(l4)

    LOGGER.info('Reads with W1 sequence found: %d / %d (%.1f %%)',
                passed_w1, c, 100*(passed_w1/float(c)))
    LOGGER.info('Read with barcodes found: %d / %d (%.1f %%)',
                passed_bc, passed_w1, 100*(passed_bc / float(passed_w1)))
    LOGGER.info('Total fraction of reads passed: %d / %d (%.1f %%)',
                passed_bc, c, 100*(passed_bc / float(c)))

    bc1_index = pd.Index(barcodes1)
    bc2_index = pd.Index(barcodes2)
    df = pd.DataFrame(np.uint64(bc_counts), index=bc1_index, columns=bc2_index)
    df.to_csv(output_count_file, sep='\t')
    
    LOGGER.info('Done.')


def get_argument_parser():

    desc = 'Process inDrop reads.'

    parser = argparse.ArgumentParser(
        description=desc, add_help=False)

    g = parser.add_argument_group('Help')
    g.add_argument('-h', '--help', action='help',
                   help='Show this help message and exit.')

    g = parser.add_argument_group('Input and output files')

    g.add_argument(
        '-rb', '--barcode-read-file', type=str, required=True,
        help='.fastq.gz file containing the reads with barcode sequences.')

    g.add_argument(
        '-rm', '--mrna-read-file', type=str, required=True,
        help='.fastq.gz file containing the reads with mRNA sequences.')

    g.add_argument(
        '-b1', '--barcode1-file', type=str, required=True,
        help='.tsv file containing first (left) barcodes.')

    g.add_argument(
        '-b2', '--barcode2-file', type=str, required=True,
        help='.tsv file containing first (left) barcodes.')

    g.add_argument(
        '-o', '--output-dir', type=str, required=True,
        help='Output directory.')

    g = parser.add_argument_group('Counting parameters')

    g.add_argument(
        '-m', '--max-reads', type=int, required=False, default=0,
        help='Maxmimum number of reads to process. If 0, '
             'process all reads. [0]'
    )

    return parser


def main(args=None):
    
    if args is None:
        parser = get_argument_parser()
        args = parser.parse_args()
    
    barcode_read_file = args.barcode_read_file
    mrna_read_file = args.mrna_read_file
    barcode1_file = args.barcode1_file
    barcode2_file = args.barcode2_file
    #output_count_file = args.output_count_file
    #output_read_file = args.output_read_file
    output_dir = args.output_dir

    max_reads = args.max_reads
    if max_reads == 0:
        max_reads=None

    count_barcodes_cython(
        barcode_read_file, mrna_read_file,
        barcode1_file, barcode2_file,
        output_dir, logger=LOGGER,
        max_reads=max_reads)


LOGGER = get_logger()

if __name__ == '__main__':
    main()
