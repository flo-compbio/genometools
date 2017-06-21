import os
import logging

from jinja2 import Environment, PackageLoader, select_autoescape

_TEMPLATE_ENV = Environment(
    loader=PackageLoader('genometools',
                         os.path.join('data', 'gcloud', 'pipelines')),
    autoescape=select_autoescape(['html', 'xml'])
)

_LOGGER = logging.getLogger(__name__)


def run_fastqc(credentials, instance_config, instance_name,
               script_dir, input_file, output_dir, self_destruct=True,
               **kwargs):
    """Run FASTQC.

    TODO: docstring"""
    template = _TEMPLATE_ENV.get_template('fastqc.sh')
    startup_script = template.render(
        script_dir=script_dir,
        input_file=input_file,
        output_dir=output_dir,
        self_destruct=self_destruct)

    if len(startup_script) > 32768:
        raise ValueError('Startup script larger than 32,768 bytes!')

    #print(startup_script)

    instance_config.create_instance(
        credentials, instance_name, startup_script=startup_script, **kwargs)


def sra_download_paired_end(credentials, instance_config, instance_name,
                            script_dir, sra_run_acc, output_dir, **kwargs):    
    """Download paired-end reads from SRA and convert to gzip'ed FASTQ files.

    TODO: docstring"""

    template = _TEMPLATE_ENV.get_template('sra_download_paired_end.sh')
    startup_script = template.render(
        script_dir=script_dir,
        sra_run_acc=sra_run_acc,
        output_dir=output_dir)

    if len(startup_script) > 32768:
        raise ValueError('Startup script larger than 32,768 bytes!')

    #print(startup_script)

    instance_config.create_instance(
        credentials, instance_name, startup_script=startup_script, **kwargs)


def trim_fastq(credentials, instance_config, instance_name,
               script_dir, input_file, output_file,
               trim_crop, trim_headcrop=0,
               **kwargs):
    """Trims a FASTQ file.
    
    TODO: docstring"""
    template = _TEMPLATE_ENV.get_template('trim_fastq.sh')
    startup_script = template.render(
        script_dir=script_dir,
        input_file=input_file,
        output_file=output_file,
        trim_crop=trim_crop,
        trim_headcrop=trim_headcrop)

    if len(startup_script) > 32768:
        raise ValueError('Startup script larger than 32,768 bytes!')

    #print(startup_script)

    instance_config.create_instance(
        credentials, instance_name, startup_script=startup_script, **kwargs)
