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

from setuptools import setup, find_packages, Command
from codecs import open
from os import path

root = 'genometools'

here = path.abspath(path.dirname(__file__))
description = 'GenomeTools: Scripts and Functions For Working With Genomic Data.'

# get long description from file
long_description = ''
with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

class CleanCommand(Command):
    """Removes files generated by setuptools.

    """
    # see https://github.com/trigger/trigger/blob/develop/setup.py
    user_options = []
    def initialize_options(self):
        pass
    def finalize_options(self):
        pass
    def run(self):
        error_msg = 'You must run this command in the package root!'
        assert os.getcwd() == here, error_msg
        os.system ('rm -rf ./dist ./build ./*.egg-info ')

setup(
    name='genometools',

    version='1.2rc1',

    description=description,
    long_description=long_description,

    # homepage
    url='https://github.com/flo-compbio/genometools',

    author='Florian Wagner',
    author_email='florian.wagner@duke.edu',

    license='GPLv3',

    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        'Development Status :: 3 - Alpha',

        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',

        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',

        'Programming Language :: Python :: 2.7',
    ],

    keywords='genome genes tools analysis',

    #packages=find_packages(exclude=['contrib', 'docs', 'tests*']),
    packages=[root],

	#libraries = [],

    install_requires = [],

    extras_require = {
            'docs': ['sphinx','sphinx-rtd-theme','sphinx-argparse','mock']
    },

	# data
    #package_data={},

	# data outside the package
    #data_files=[('my_data', ['data/data_file'])],

	# executable scripts
    entry_points = {
        'console_scripts': [
            'extract_protein_coding_genes.py = genometools.extract_protein_coding_genes:main',
            'extract_entrez2gene.py = genometools.extract_entrez2gene:main',
            'sra_download_experiment.py = genometools.sra.download_experiment:main',
        ],
    },

    cmdclass = {
        'clean': CleanCommand,
    },

)
