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

def get_logger(log_file=None,log_level=logging.INFO):
    """Returns a fully configured `logging.Logger` object.

    This function serves to simplify the generation of a logger that either
    writes either to a file or to stdout.

    Parameters
    ----------
    log_file: Optional[str]
        The path of the log file. If None, write to stdout.
    log_level: Optional[int]
        A logging level as 
        `defined <https://docs.python.org/2/library/logging.html#logging-levels>`_ in
        Python's logging module.

    Returns
    -------
    logger: logging.Logger
        The logger object.

    """

    log_format = '[%(asctime)s] %(levelname)s: %(message)s'
    log_datefmt = '%Y-%m-%d %H:%M:%S'
    logging.basicConfig(filename=log_file,stream=sys.stdout,level=log_level,format=log_format,datefmt=log_datefmt)
    logger = logging.getLogger()
    return logger

def open_plain_or_gzip(fn):
    """Returns a handle for a file that is either gzip'ed or not.

    Parameters
    ----------
    fn: str
        The path of the file.

    Returns
    -------
    file
        A file object.

    Notes
    -----
    Generally, reading gzip with gzip.open is very slow, and it is preferable
    to pipe the file into the python script (so that the scripts reads the file
    from stdin).
    """

    try:
        gzip.open(fn,'rb').next()
    except IOError:
        return open(fn,'r')
    else:
        return gzip.open(fn,mode)

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


