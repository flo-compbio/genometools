import argparse

from genometools import ensembl

def get_gtf_argument_parser(description):
    """Return an argument parser with basic options for reading GTF files.

    Parameters
    ----------
    description: str
        Description of the ArgumentParser

    Returns
    -------
    `argparse.ArgumentParser` object
        The argument parser.
    """
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('-a','--annotation-file', default='-',
        help="""Path of Ensembl gene annotation file (in GTF format). The file
                may be gzip'ed. If set to ``-``, read from ``stdin``.""")

    parser.add_argument('-o','--output-file', required=True,
        help="""Path of output file. If set to ``-``, print to ``stdout``,
                and redirect logging messages to ``stderr``.""")

    parser.add_argument('-s', '--species',
        choices=sorted(ensembl.species_chrompat.keys()), default='human',
        help="""Species for which to extract genes. (This parameter is ignored
                if ``--chromosome-pattern`` is specified.)""")

    parser.add_argument('-c', '--chromosome-pattern',
        required=False, default=None,
        help="""Regular expression that chromosome names have to match.
                If not specified, determine pattern based on ``--species``.""")

    parser.add_argument('-f','--field-name', default='gene',
        help="""Rows in the GTF file that do not contain this value
                in the third column are ignored.""")

    parser.add_argument('-l','--log-file', default=None,
        help='Path of log file. If not specified, print to stdout.')

    parser.add_argument('-q','--quiet', action='store_true',
        help='Suppress all output except warnings and errors.')

    parser.add_argument('-v','--verbose', action='store_true',
        help='Enable verbose output. Ignored if ``--quiet`` is specified.')

    return parser
