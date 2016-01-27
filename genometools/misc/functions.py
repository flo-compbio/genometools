# Copyright (c) 2015, 2016 Florian Wagner
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
import io
import errno
import shutil
import urllib2
import sys
import bisect
import gzip
import logging
import contextlib

import unicodecsv as csv

logger = logging.getLogger(__name__)

def try_open_gzip(path):
    fh = None
    try:
        gzip.open(path).next()
    except IOError:
        pass
    else:
        fh = gzip.open(path)
    return fh

@contextlib.contextmanager
def smart_open_read(path = None, mode = 'r', encoding = None, try_gzip = False):
    """Open a file for reading or return ``stdin``.

    Adapted from StackOverflow user "Wolph"
    (http://stackoverflow.com/a/17603000).
    """

    assert mode in ('r', 'rb')
    assert path is None or isinstance(path, (str, unicode))
    assert isinstance (mode, str)
    assert encoding is None or isinstance(encoding, str)
    assert isinstance(try_gzip, bool)

    fh = None
    binfh = None
    gzfh = None
    if path is None:
        # open stdin
        fh = io.open(sys.stdin.fileno(), mode = mode, encoding = encoding)

    else:
        # open an actual file

        if try_gzip:
            # gzip.open defaults to mode 'rb'
            gzfh = try_open_gzip(path)

        if gzfh is not None:
            logger.debug('Opening gzip''ed file.')
            # wrap gzip stream
            binfh = io.BufferedReader(gzfh)
            if 'b' not in mode:
                # add a text wrapper on top
                logger.debug('Adding text wrapper.')
                fh = io.TextIOWrapper(binfh, encoding = encoding)

        else:
            fh = io.open(path, mode = mode, encoding = encoding)

    yield_fh = fh
    if fh is None:
        yield_fh = binfh

    try:
        yield yield_fh

    finally:
        # close all open files
        if fh is not None:
            # make sure we don't close stdin
            if fh.fileno() != sys.stdin.fileno():
                fh.close()

        if binfh is not None:
            binfh.close()

        if gzfh is not None:
            gzfh.close()
        
@contextlib.contextmanager
def smart_open_write(path = None, mode = 'wb', encoding = None):
    """Open a file for writing or return ``stdout``.

    Adapted from StackOverflow user "Wolph"
    (http://stackoverflow.com/a/17603000).
    """
    if path is not None:
        # open a file
        fh = io.open(path, mode = mode, encoding = encoding)
    else:
        # open stdout
        fh = io.open(sys.stdout.fileno(), mode = mode, encoding = encoding)
        #fh = sys.stdout

    try:
        yield fh

    finally:
        # make sure we don't close stdout
        if fh.fileno() != sys.stdout.fileno():
            fh.close()

def test_dir_writable(path):
    """Test if we can write to a directory.
    """
    dir_ = os.path.dirname(path)
    if dir_ == '':
        dir_ = '.'
    return os.access(dir_, os.W_OK)

def test_file_writable(path):
    """Test if we can write to a file."""
    if os.path.isfile(path):
        # file exists, can we modify it?
        try:
            with open(path, 'ab') as fh:
                pass
        except IOError:
            return False
        else:
            return True
    else:
        # file does not exist, can we write to the directory?
        return test_dir_writable(path)

def download_url(url, download_file):
    """Downloads a file from a given URL.

    Source: StackOverflow user "J.F. Sebastian"
    (http://stackoverflow.com/a/11768443)

    Parameters
    ----------
    url: str
        The URL (source).
    output_file: str
        The path of the output file.
    """
    with contextlib.closing(urllib2.urlopen(url)) as fh:
        with open(download_file, 'wb') as ofh:
            shutil.copyfileobj(fh, ofh)

def make_sure_dir_exists(d, create_subfolders = False):
    """Ensures that a directory exists.

    Adapted from StackOverflow users "Bengt" and "Heikki Toivonen"
    (http://stackoverflow.com/a/5032238).
    """
    try:
        if create_subfolders:
            os.makedirs(d)
        else:
            os.mkdir(d)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

def open_plain_or_gzip(fn, mode = 'rb'):
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
        gzip.open(fn, 'rb').next()
    except IOError:
        return io.open(fn, mode)
    else:
        return gzip.open(fn, 'rb')

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

def read_single(fn, enc = 'utf-8'):
    """ Reads the first column of a tab-delimited text file.

    The file can either be uncompressed or gzip'ed.

    Parameters
    ----------
    fn: str
        The path of the file.
    enc: str, optional
        The file encoding.

    Returns
    -------
    List[str]
        A list containing the elements in the first column.

    """
    data = []
    with open_plain_or_gzip(fn) as fh:
        reader = csv.reader(fh, dialect = 'excel-tab', encoding = enc)
        for l in reader:
            data.append(l[0])
    return data

def read_all(fn, enc = 'utf-8'):
    """ Reads a tab-delimited text file.

    The file can either be uncompressed or gzip'ed.

    Parameters
    ----------
    fn: str
        The path of the file.
    enc: str, optional
        The file encoding.

    Returns
    -------
    List[List[str]]
        A list, which each element containing the contents of a row
        (as a list).
    """
    data = []
    with open_plain_or_gzip(fn) as fh:
        reader = csv.reader(fh, dialect='excel-tab', encoding = enc)
        for l in reader:
            data.append(l)
    return data
