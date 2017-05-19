import os
import logging

from jinja2 import Environment, PackageLoader, select_autoescape

_TEMPLATE_ENV = Environment(
    loader=PackageLoader('genometools',
                         os.path.join('data', 'gcloud', 'pipelines',
                                      'inDrop')),
    autoescape=select_autoescape(['html', 'xml'])
)

_LOGGER = logging.getLogger(__name__)


def process_reads(credentials, instance_config, instance_name,
                  script_dir,
                  barcode_read_file, mrna_read_file,
                  barcode1_file, barcode2_file,
                  output_dir,
                  max_reads=None,
                  self_destruct=True,
                  **kwargs):
    """Process inDrop reads.

    Recommended machine type: "n1-standard-2" (7.5 GB of RAM, 2 vCPUs)

    TODO: docstring"""

    if max_reads is None:
        max_reads = 0

    template = _TEMPLATE_ENV.get_template(
        os.path.join('process_reads.sh'))
    startup_script = template.render(
        script_dir=script_dir,
        barcode_read_file=barcode_read_file,
        mrna_read_file=mrna_read_file,
        barcode1_file=barcode1_file,
        barcode2_file=barcode2_file,
        output_dir=output_dir,
        max_reads=max_reads,
        self_destruct=self_destruct)

    if len(startup_script) > 32768:
        raise ValueError('Startup script larger than 32,768 bytes!')

    #print(startup_script)

    op_name = instance_config.create_instance(
        credentials, instance_name, startup_script=startup_script, **kwargs)

    return op_name


def run_pipeline(credentials, instance_config, instance_name,
                 script_dir,
                 barcode_read_file, mrna_read_file,
                 barcode1_file, barcode2_file,
                 index_dir,
                 output_dir,
                 max_reads=None,
                 self_destruct=True,
                 num_threads=16,
                 **kwargs):

    """Process inDrop reads.

    Recommended machine type: "n1-standard-2" (7.5 GB of RAM, 2 vCPUs)

    TODO: docstring"""

    if max_reads is None:
        max_reads = 0

    template = _TEMPLATE_ENV.get_template(
        os.path.join('indrop.sh'))
    startup_script = template.render(
        script_dir=script_dir,
        barcode_read_file=barcode_read_file,
        mrna_read_file=mrna_read_file,
        barcode1_file=barcode1_file,
        barcode2_file=barcode2_file,
        index_dir=index_dir,
        output_dir=output_dir,
        max_reads=max_reads,
        num_threads=num_threads,
        self_destruct=self_destruct)

    if len(startup_script) > 32768:
        raise ValueError('Startup script larger than 32,768 bytes!')

    #print(startup_script)

    op_name = instance_config.create_instance(
        credentials, instance_name, startup_script=startup_script, **kwargs)

    return op_name    
    