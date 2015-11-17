import pkg_resources

__version__ = pkg_resources.require('genometools')[0].version

__all__ = ['misc']
