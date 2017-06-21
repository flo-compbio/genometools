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

"""Miscellaneous functions that are useful in many different contexts.

"""

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
_oldstr = str
from builtins import *

import os
import io
import errno
import shutil
import sys
import bisect
import gzip
import logging
import contextlib
import subprocess as subproc
import locale
import ftplib
import hashlib

import six

if six.PY3:
    from urllib import parse as urlparse
else:
    import urlparse


import unicodecsv as csv
import requests

logger = logging.getLogger(__name__)


def try_open_gzip(path):
    fh = None
    try:
        next(gzip.open(path))
    except IOError:
        pass
    else:
        fh = gzip.open(path)
    return fh


def get_file_md5sum(path):
    """Calculate the MD5 hash for a file."""
    with open(path, 'rb') as fh:
        h = str(hashlib.md5(fh.read()).hexdigest())
    return h


@contextlib.contextmanager
def smart_open_read(path=None, mode='rb', encoding=None, try_gzip=False):
    """Open a file for reading or return ``stdin``.

    Adapted from StackOverflow user "Wolph"
    (http://stackoverflow.com/a/17603000).
    """
    assert mode in ('r', 'rb')
    assert path is None or isinstance(path, (str, _oldstr))
    assert isinstance(mode, (str, _oldstr))
    assert encoding is None or isinstance(encoding, (str, _oldstr))
    assert isinstance(try_gzip, bool)

    fh = None
    binfh = None
    gzfh = None
    if path is None:
        # open stdin
        fh = io.open(sys.stdin.fileno(), mode=mode, encoding=encoding)

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
                fh = io.TextIOWrapper(binfh, encoding=encoding)

        else:
            fh = io.open(path, mode=mode, encoding=encoding)

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
def smart_open_write(path=None, mode='wb', encoding=None):
    """Open a file for writing or return ``stdout``.

    Adapted from StackOverflow user "Wolph"
    (http://stackoverflow.com/a/17603000).
    """
    if path is not None:
        # open a file
        fh = io.open(path, mode=mode, encoding=encoding)
    else:
        # open stdout
        fh = io.open(sys.stdout.fileno(), mode=mode, encoding=encoding)
        #fh = sys.stdout

    try:
        yield fh

    finally:
        # make sure we don't close stdout
        if fh.fileno() != sys.stdout.fileno():
            fh.close()

def test_dir_writable(path):
    """Test if we can write to a directory.

    Parameters
    ----------
    dir: str
        The directory path.

    Returns
    -------
    bool
        Whether the directory is writable or not.
    """
    dir_ = os.path.dirname(path)
    if dir_ == '':
        dir_ = '.'
    return os.access(dir_, os.W_OK)

def test_file_writable(path):
    """Test if we can write to a file.

    Parameters
    ----------
    path: str
        The file path.

    Returns
    -------
    bool
        Whether the file is writable or not.
    """
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

def get_fize_size(path):
    """The the size of a file.

    Parameters
    ----------
    path: str
        The file path.

    Returns
    -------
    int
        The size of the file in bytes.
    """
    return os.path.getsize(path)

def get_url_size(url):
    """Get the size of a URL.

    Note: Uses requests, so it does not work for FTP URLs.

    Source: StackOverflow user "Burhan Khalid".
    (http://stackoverflow.com/a/24585314/5651021)

    Parameters
    ----------
    url : str
        The URL.

    Returns
    -------
    int
        The size of the URL in bytes.
    """
    r = requests.head(url, headers={'Accept-Encoding': 'identity'})
    size = int(r.headers['content-length'])
    return size


def get_url_file_name(url):
    """Get the file name from an url
    
    Parameters
    ----------
    url : str

    Returns
    -------
    str
        The file name 
    """

    assert isinstance(url, (str, _oldstr))
    return urlparse.urlparse(url).path.split('/')[-1]

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
    assert isinstance(dir_, (str, _oldstr))
    assert isinstance(create_subfolders, bool)

    try:
        if create_subfolders:
            os.makedirs(dir_)
        else:
            os.mkdir(dir_)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

def get_file_size(path):
    """The the size of a file in bytes.

    Parameters
    ----------
    path: str
        The path of the file.

    Returns
    -------
    int
        The size of the file in bytes.

    Raises
    ------
    IOError
        If the file does not exist.
    OSError
        If a file system error occurs.
    """
    assert isinstance(path, (str, _oldstr))

    if not os.path.isfile(path):
        raise IOError('File "%s" does not exist.', path)

    return os.path.getsize(path)

def get_file_checksum(path):
    """Get the checksum of a file (using ``sum``, Unix-only).

    This function is only available on certain platforms.

    Parameters
    ----------
    path: str
        The path of the file.

    Returns
    -------
    int
        The checksum.

    Raises
    ------
    IOError
        If the file does not exist.
    """

    if not (sys.platform.startswith('linux') or \
                        sys.platform in ['darwin', 'cygwin']):
        raise OSError('This function is not available on your platform.')

    assert isinstance(path, (str, _oldstr))

    if not os.path.isfile(path): # not a file
        raise IOError('File "%s" does not exist.' %(path))

    # calculate checksum
    sub = subproc.Popen('sum "%s"' %(path), bufsize=-1, shell=True,
                        stdout=subproc.PIPE)
    stdoutdata = sub.communicate()[0]
    assert sub.returncode == 0

    # in Python 3, communicate() returns bytes that need to be decoded
    encoding = locale.getpreferredencoding()
    stdoutstr = str(stdoutdata, encoding=encoding)

    file_checksum = int(stdoutstr.split(' ')[0])
    logger.debug('Checksum of file "%s": %d', path, file_checksum)
    return file_checksum


def test_file_checksum(path, checksum):
    """Test if a file has a given checksum (using ``sum``, Unix-only).

    Parameters
    ----------
    path: str
        The path of the file.
    checksum: int
        The checksum to compare.
    
    Returns
    -------
    bool
        Whether or not the file has the given checksum.

    Raises
    ------
    IOError
        If the file does not exist.
    """
    assert isinstance(path, (str, _oldstr))
    assert isinstance(checksum, int)

    # calculate file checksum and compare to given checksum
    file_checksum = get_file_checksum(path)
    logger.debug('File checksum: %d. Reference checksum: %d. Match: %s.',
                 file_checksum, checksum, str(file_checksum == checksum))
    return file_checksum == checksum


# @contextlib.contextmanager
def gzip_open_text(path, encoding=None):
    """Opens a plain-text file that may be gzip'ed.

    Parameters
    ----------
    path : str
        The file.
    encoding : str, optional
        The encoding to use.

    Returns
    -------
    file-like
        A file-like object.

    Notes
    -----
    Generally, reading gzip'ed files with gzip.open is very slow, and it is
    preferable to pipe the file into the python script using ``gunzip -c``.
    The script then reads the file from stdin.
    """
    if encoding is None:
        encoding = sys.getdefaultencoding()

    assert os.path.isfile(path)

    is_compressed = False
    try:
        gzip.open(path, mode='rb').read(1)
    except IOError:
        pass
    else:
        is_compressed = True

    if is_compressed:
        if six.PY2:
            import codecs
            zf = gzip.open(path, 'rb')
            reader = codecs.getreader(encoding)
            fh = reader(zf)

        else:
            fh = gzip.open(path, mode='rt', encoding=encoding)

    else:
        # the following works in Python 2.7, thanks to future
        fh = open(path, mode='r', encoding=encoding)

    return fh


def is_writable(path):
    """Tests if a file is writable."""
    try:
        with open(path, 'a'):
            pass
    except:
        return False
    return True


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


def read_single(path, encoding = 'UTF-8'):
    """ Reads the first column of a tab-delimited text file.

    The file can either be uncompressed or gzip'ed.

    Parameters
    ----------
    path: str
        The path of the file.
    enc: str
        The file encoding.

    Returns
    -------
    List of str
        A list containing the elements in the first column.

    """
    assert isinstance(path, (str, _oldstr))
    data = []
    with smart_open_read(path, mode='rb', try_gzip=True) as fh:
        reader = csv.reader(fh, dialect='excel-tab', encoding=encoding)
        for l in reader:
            data.append(l[0])
    return data


def read_all(path, encoding='UTF-8'):
    """ Reads a tab-delimited text file.

    The file can either be uncompressed or gzip'ed.

    Parameters
    ----------
    path: str
        The path of the file.
    enc: str, optional
        The file encoding.

    Returns
    -------
    List of (tuple of str)
        A list, which each element containing the contents of a row
        (as a tuple).
    """
    assert isinstance(path, (str, _oldstr))
    data = []
    with smart_open_read(path, mode='rb', try_gzip=True) as fh:
        reader = csv.reader(fh, dialect='excel-tab', encoding=encoding)
        for l in reader:
            data.append(tuple(l))
    return data
