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

"""Miscellaneous functions that are useful in many different contexts.

"""

import os
import sys
import csv
import bisect
import gzip
import logging
import contextlib

def configure_logger(name, log_stream = sys.stdout, log_file = None,
        log_level = logging.INFO, keep_old_handlers = False,
        propagate = False):
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

    __ log_lvl_

    Returns
    -------
    logger: `logging.Logger`
        The logger.

    Notes
    -----
    Note that if ``log_stream`` and ``log_file`` are both ``None``, no handlers
    will be created.

    .. _log_lvl: https://docs.python.org/2/library/logging.html#logging-levels

    """

    # create a child logger
    logger = logging.getLogger(name)

    # set the logger's level
    logger.setLevel(log_level)

    # set the logger's propagation attribute
    logger.propgate = propagate

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
      
    return logger

@contextlib.contextmanager
def smart_open(filename=None,mode='r',try_gzip=False):
    """Open a file for reading or return ``stdin``.

    Author: StackOverflow user "Wolph"
    Source: http://stackoverflow.com/a/17603000
    """
    fh = None
    if filename and filename != '-':
        if try_gzip:
            fh = open_plain_or_gzip(filename,mode)
        else:
            fh = open(filename, mode)
    else:
        fh = sys.stdin

    try:
        yield fh
    finally:
        if fh is not sys.stdin:
            fh.close()

@contextlib.contextmanager
def smart_open_write(filename=None,mode='w'):
    """Open a file for writing or return ``stdout``.

    Author: StackOverflow user "Wolph"
    Source: http://stackoverflow.com/a/17603000
    """
    if filename and filename != '-':
        fh = open(filename, mode)
    else:
        fh = sys.stdout

    try:
        yield fh
    finally:
        if fh is not sys.stdout:
            fh.close()

def open_plain_or_gzip(fn,mode='r'):
    """Returns a handle for a file that is either gzip'ed or not.

    Parameters
    ----------
    fn: str
        The path of the file.
    mode: str
        The mode to be used in opening the file if it is not gzip'ed.

    Returns
    -------
    file
        A file object.

    Notes
    -----
    Generally, reading gzip'ed files with gzip.open is very slow, and it is
    preferable to pipe the file into the python script using ``gunzip -c``.
    The script then reads the file from stdin.

    """

    try:
        gzip.open(fn,'rb').next()
    except IOError:
        return open(fn,mode)
    else:
        return gzip.open(fn,'rb')

def flatten(l):
    """Flattens a list of lists.

    Parameters
    ----------
    l: list
        The list of lists.

    Returns
    -------
    list
        The flattened list.

    """
    # see http://stackoverflow.com/questions/952914/making-a-flat-list-out-of-list-of-lists-in-python#comment10547502_952952
    # use incomprensible list comprehension
    return [item for sublist in l for item in sublist]

def bisect_index(a, x):
    """ Find the leftmost index of an element in a list using binary search.

    Parameters
    ----------
    a: list
        A sorted list.
    x: arbitrary
        The element.

    Returns
    -------
    int
        The index.

    """

    i = bisect.bisect_left(a, x)
    if i != len(a) and a[i] == x:
        return i
    raise ValueError

def argsort(seq):
    """ Returns a list of indices that would sort a list.

    Parameters
    ----------
    seq: List
        The list.

    Returns
    -------
    List[int]
        The list of indices that would sort the given list ``seq``.

    Notes
    -----
    If the returned list of indices can be a NumPy array, use `numpy.lexsort`
    instead. If the given list ``seq`` is a NumPy array, use `numpy.argsort`
    instead.

    """

    # see http://stackoverflow.com/questions/3382352/equivalent-of-numpy-argsort-in-basic-python/3382369#3382369
    return sorted(range(len(seq)), key=seq.__getitem__)

def argmin(seq):
    """ Obtains the index of the smallest element in a list.

    Parameters
    ----------
    seq: List
        The list.

    Returns
    -------
    int
        The index of the smallest element.

    """

    return argsort(seq)[0]

def argmax(seq):
    """ Obtains the index of the largest element in a list.

    Parameters
    ----------
    seq: List
        The list

    Returns
    -------
    int
        The index of the largest element.
    """

    return argsort(seq)[-1]

def read_single(fn):
    """ Reads the first column of a tab-delimited text file.

    The file can either be uncompressed or gzip'ed.

    Parameters
    ----------
    fn: str
        The path of the file.

    Returns
    -------
    List[str]
        A list containing the elements in the first column.

    """
    data = []
    with open_plain_or_gzip(fn) as fh:
        reader = csv.reader(fh,dialect='excel-tab')
        for l in reader:
            data.append(l[0])
    return data

def read_all(fn):
    """ Reads a tab-delimited text file.

    The file can either be uncompressed or gzip'ed.

    Parameters
    ----------
    fn: str
        The path of the file.

    Returns
    -------
    List[List[str]]
        A list, which each element containing the contents of a row
        (as a list).
    """

    data = []
    with open_plain_or_gzip(fn) as fh:
        reader = csv.reader(fh,dialect='excel-tab')
        for l in reader:
            data.append(l)
    return data


