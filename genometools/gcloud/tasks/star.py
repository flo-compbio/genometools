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
                   num_threads=8, seed_start_lmax=50,
                   mismatch_nmax=10, multimap_nmax=10, splice_min_overhang=3,
                   out_mult_nmax=-1, sort_bam=True,
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

    template = _TEMPLATE_ENV.get_template(
        os.path.join('map_single-end.sh'))
    startup_script = template.render(
        script_dir=script_dir,
        index_dir=index_dir,
        fastq_file=fastq_file,
        output_dir=output_dir,
        num_threads=num_threads,
        seed_start_lmax=seed_start_lmax,
        self_destruct=self_destruct,
        mismatch_nmax=mismatch_nmax,
        multimap_nmax=multimap_nmax,
        splice_min_overhang=splice_min_overhang,
        out_mult_nmax=out_mult_nmax,
        compressed=compressed,
        out_sam_type=out_sam_type)

    if len(startup_script) > 32768:
        raise ValueError('Startup script larger than 32,768 bytes!')

    #print(startup_script)

    op_name = instance_config.create_instance(
        credentials, instance_name, startup_script=startup_script, **kwargs)

    return op_name


def map_paired_end(credentials, instance_config, instance_name,
                   script_dir, index_dir, fastq_file1, fastq_file2, output_dir,
                   num_threads=1, max_mate_dist=1000000, chim_segment_min=0,
                   seed_start_lmax=50,
                   score_min_abs=0, score_min_rel=0.66,
                   match_nmin_abs=0, match_nmin_rel=0.66,
                   **kwargs):
    """Maps single-end reads using STAR.

    Recommended machine type: "n1-highmem-8" (52 GB of RAM, 8 vCPUs)

    TODO: docstring"""

    template = _TEMPLATE_ENV.get_template(
        os.path.join('map_paired-end.sh'))
    startup_script = template.render(
        script_dir=script_dir,
        index_dir=index_dir,
        fastq_file1=fastq_file1,
        fastq_file2=fastq_file2,
        output_dir=output_dir,
        num_threads=num_threads,
        max_mate_dist=max_mate_dist,
        chim_segment_min=chim_segment_min,
        seed_search_lmax=seed_start_lmax,
        score_min_abs=score_min_abs,
        score_min_rel=score_min_rel,
        match_nmin_abs=match_nmin_abs,
        match_nmin_rel=match_nmin_rel)

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
