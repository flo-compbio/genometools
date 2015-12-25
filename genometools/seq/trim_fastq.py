#!/usr/bin/env python2.7

import sys
import argparse
import itertools as it

from genometools import cli

def get_argument_parser():

    desc = 'Trim FASTQ file (read from stdin, write to stdout).'
    parser = cli.get_argument_parser(desc = desc)

    parser.add_argument('-l','--left', type=int, default=0,
        help="""Number of base pairs to trim from the left.""")

    parser.add_argument('-r','--right', type=int, default=0,
        help="""Number of base pairs to trim from the right.""")

    return parser

def main(args=None):

    if args is None:
        parser = get_argument_parser()
        args = parser.parse_args()

    left = args.left
    right = args.right

    #counter = it.cycle(range(4))
    from_right = right + 1
    while True:
        try:
            sys.stdout.write(sys.stdin.next())
        except StopIteration:
            break
        sys.stdout.write(sys.stdin.next()[left:-from_right] + '\n')
        sys.stdout.write(sys.stdin.next())
        sys.stdout.write(sys.stdin.next()[left:-from_right] + '\n')

    return 0

if __name__ == '__main__':
    return_code = main()
    sys.exit(return_code)
