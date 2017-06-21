# Copyright (c) 2017 Florian Wagner
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

"""Functions for running STAR on Google Cloud."""

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
_oldstr = str
from builtins import *
# import six

import os
import logging

from jinja2 import Environment, PackageLoader, select_autoescape

_TEMPLATE_ENV = Environment(
    loader=PackageLoader('genometools',
                         os.path.join('data', 'gcloud', 'pipelines', 'STAR')),
    autoescape=select_autoescape(['html', 'xml'])
)

_LOGGER = logging.getLogger(__name__)


def map_single_end(credentials, instance_config, instance_name,
                   script_dir, index_dir, fastq_file, output_dir,
                   num_threads=None, seed_start_lmax=None,
                   mismatch_nmax=None, multimap_nmax=None,
                   splice_min_overhang=None,
                   out_mult_nmax=None, sort_bam=True, keep_unmapped=False,
                   self_destruct=True, compressed=True,
                   **kwargs):
    """Maps single-end reads using STAR.

    Reads are expected in FASTQ format. By default, they are also expected to
    be compressed with gzip.

    - recommended machine type: "n1-standard-16" (60 GB of RAM, 16 vCPUs).
    - recommended disk size: depends on size of FASTQ files, at least 128 GB.

    TODO: docstring"""

    if sort_bam:
        out_sam_type = 'BAM SortedByCoordinate'
    else:
        out_sam_type = 'BAM Unsorted'

    # template expects a list of FASTQ files
    fastq_files = fastq_file
    if isinstance(fastq_files, (str, _oldstr)):
        fastq_files = [fastq_file]

    template = _TEMPLATE_ENV.get_template(
        os.path.join('map_single-end.sh'))
    startup_script = template.render(
        script_dir=script_dir,
        index_dir=index_dir,
        fastq_files=fastq_files,
        output_dir=output_dir,
        num_threads=num_threads,
        seed_start_lmax=seed_start_lmax,
        self_destruct=self_destruct,
        mismatch_nmax=mismatch_nmax,
        multimap_nmax=multimap_nmax,
        splice_min_overhang=splice_min_overhang,
        out_mult_nmax=out_mult_nmax,
        keep_unmapped=keep_unmapped,
        compressed=compressed,
        out_sam_type=out_sam_type)

    if len(startup_script) > 32768:
        raise ValueError('Startup script larger than 32,768 bytes!')

    #print(startup_script)

    op_name = instance_config.create_instance(
        credentials, instance_name, startup_script=startup_script, **kwargs)

    return op_name


def map_paired_end(credentials, instance_config, instance_name,
                   script_dir, index_dir, fastq_file_pair, output_dir,
                   num_threads=None, 
                   max_intron_size=None, max_mate_dist=None,
                   splice_min_overhang=None,
                   out_mult_nmax=None, sort_bam=True, keep_unmapped=False,
                   seed_start_lmax=None,
                   mismatch_nmax=None, multimap_nmax=None,
                   self_destruct=True, compressed=True,
                   **kwargs):
                   #seed_start_lmax=50,
                   #score_min_abs=0, score_min_rel=0.66,
                   #match_nmin_abs=0, match_nmin_rel=0.66,
    """Maps single-end reads using STAR.

    Recommended machine type: "n1-highmem-8" (52 GB of RAM, 8 vCPUs)

    TODO: docstring"""

    if isinstance(fastq_file_pair[0], (str, _oldstr)):
        fastq_file_pairs = [fastq_file_pair]
    else:
        fastq_file_pairs = fastq_file_pair

    if sort_bam:
        out_sam_type = 'BAM SortedByCoordinate'
    else:
        out_sam_type = 'BAM Unsorted'

    template = _TEMPLATE_ENV.get_template(
        os.path.join('map_paired-end.sh'))
    startup_script = template.render(
        script_dir=script_dir,
        index_dir=index_dir,
        fastq_file_pairs=fastq_file_pairs,
        output_dir=output_dir,
        compressed=compressed,
        num_threads=num_threads,
        seed_start_lmax=seed_start_lmax,
        self_destruct=self_destruct,
        mismatch_nmax=mismatch_nmax,
        multimap_nmax=multimap_nmax,
        splice_min_overhang=splice_min_overhang,
        out_mult_nmax=out_mult_nmax,
        out_sam_type=out_sam_type,
        max_intron_size=max_intron_size,
        max_mate_dist=max_mate_dist,
        keep_unmapped=keep_unmapped)
        #chim_segment_min=chim_segment_min,
        #seed_search_lmax=seed_start_lmax,
        #score_min_abs=score_min_abs,
        #score_min_rel=score_min_rel,
        #match_nmin_abs=match_nmin_abs,
        #match_nmin_rel=match_nmin_rel)

    if len(startup_script) > 32768:
        raise ValueError('Startup script larger than 32,768 bytes!')

    #print(startup_script)

    op_name = instance_config.create_instance(
        credentials, instance_name, startup_script=startup_script, **kwargs)

    return op_name


def generate_index(credentials, instance_config, instance_name,
                   script_dir, genome_file, output_dir, annotation_file=None,
                   splice_overhang=100,
                   num_threads=8, chromosome_bin_bits=18,
                   genome_memory_limit=31000000000,
                   self_destruct=True,
                   **kwargs):
    """Generates a STAR index.

    Recommended machine type: "n1-highmem-8" (52 GB of RAM, 8 vCPUs)

    TODO: docstring"""

    template = _TEMPLATE_ENV.get_template(
        os.path.join('generate_index.sh'))
    startup_script = template.render(
        script_dir=script_dir,
        genome_file=genome_file,
        annotation_file=annotation_file,
        splice_overhang=splice_overhang,
        output_dir=output_dir,
        num_threads=num_threads,
        chromosome_bin_bits=chromosome_bin_bits,
        genome_memory_limit=genome_memory_limit,
        self_destruct=self_destruct)

    op_name = instance_config.create_instance(
        credentials, instance_name, startup_script=startup_script, **kwargs)

    return op_name
