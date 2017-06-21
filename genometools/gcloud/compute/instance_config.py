"""Module containing the `InstanceConfig` class."""

from . import create_instance

class InstanceConfig(object):
    """Class that holds instance configuration data.
    
    TODO: docstring"""

    def __init__(self, project, zone,
                 machine_type='f1-micro', disk_size_gb=10):
        self.project = project
        self.zone = zone
        self.machine_type = machine_type
        self.disk_size_gb = disk_size_gb


    def create_instance(self, credentials, name, **kwargs):
        """Create an instance based on the configuration data.
        
        TODO: docstring"""
        op_name = create_instance(
            credentials, self.project, self.zone, name,
            machine_type=self.machine_type,
            disk_size_gb=self.disk_size_gb,
            **kwargs)
        
        return op_name
