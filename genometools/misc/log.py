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
import logging

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

def add_logging_params(parser):
    """Adds logging parameters to an argument parser.

    Parameters
    ----------
    parser: `argparse.ArgumentParser`
        The argument parser.

    Returns:
    `argparse.ArgumentParser`
        The argument parser with the parameters added.
    """

    rep = parser.add_argument_group('Reporting options')
    rep.add_argument('-l','--log-file', default=None,
        help='Path of log file. If not specified, print to stdout.')

    rep.add_argument('-q','--quiet', action='store_true',
        help='Suppress all output except warnings and errors.')

    rep.add_argument('-v','--verbose', action='store_true',
        help='Enable verbose output. Ignored if ``--quiet`` is specified.')

    return parser
