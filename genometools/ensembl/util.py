import textwrap

from genometools import cli
from genometools import ensembl

def get_gtf_argument_parser(desc, default_field_name = 'gene'):
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
    parser = cli.get_argument_parser(desc = desc)

    str_type = cli.str_type

    parser.add_argument('-a','--annotation-file', default='-', type = str_type,
            help = textwrap.dedent("""\
                Path of Ensembl gene annotation file (in GTF format). The file
                may be gzip'ed. If set to ``-``, read from ``stdin``."""))

    parser.add_argument('-o','--output-file', required=True, type = str_type,
            help = textwrap.dedent("""\
                Path of output file. If set to ``-``, print to ``stdout``,
                and redirect logging messages to ``stderr``."""))

    parser.add_argument('-s', '--species', type = str_type,
            choices = sorted(ensembl.species_chrompat.keys()), default = 'human',
            help = textwrap.dedent("""\
                Species for which to extract genes. (This parameter is ignored
                if ``--chromosome-pattern`` is specified.)"""))

    parser.add_argument('-c', '--chromosome-pattern', type = str_type,
            required = False, default = None, help = textwrap.dedent("""\
                Regular expression that chromosome names have to match.
                If not specified, determine pattern based on
                ``--species``."""))

    parser.add_argument('-f', '--field-name', default = default_field_name,
            type = str_type, help = textwrap.dedent("""\
                Rows in the GTF file that do not contain this value
                in the third column are ignored."""))

    cli.add_reporting_args(parser)

    return parser
