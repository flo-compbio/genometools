import os
import re

from genometools import misc

class EnsemblGeneAnnotationData(object):
    
    file_pat = re.compile(r'(.*)?\.(.*)\.(\d+)\.gtf\.gz$')
    
    def __init__(self, dir_, species, middle, release):
        self.dir_ = dir_.rstrip(os.sep)
        self.species = species
        self.middle = middle
        self.release = release
    
    def __repr__(self):
        return '<EnsemblGeneAnnotationData (dir_=%s; species=%s; middle=%s; release=%d)>' \
                %(self.species, self.middle, self.release)

    def __str__(self):
        return '<EnsemblGeneAnnotationData for species "%s", release %d (dir_="%s")>' \
                %(self.species, self.release, self.dir_)

    @classmethod
    def from_path(cls, path):
        dir_, file_name = os.path.split(path)
        m = cls.file_pat.match(file_name)
        species = m.group(1)
        middle = m.group(2)
        release = int(m.group(3))
        return cls(dir_, species, middle, release)

    @property
    def file_name(self):
        return '%s.%s.%d.gtf.gz' %(self.species, self.middle, self.release)
    
    @property
    def path(self):
        return '%s/%s' %(self.dir_, self.file_name)

    @classmethod
    def scan_dir(cls, dir_):
        dir_ = dir_.rstrip(os.sep)
        data = []
        for f in os.listdir(dir_):
            path = dir_ + os.sep + f
            if os.path.isfile(path) and cls.file_pat.match(f) is not None:
                data.append(cls.from_path(path))

        data = sorted(data, key = lambda d: [d.species, d.release])
        return data
    
class UniProtGOAnnotationData(object):
    
    file_pat = re.compile(r'gene_association.goa_(.*?).(%d+).gz')
    ontology_file_pat = re.compile(r'go-basic_(\d{4}-\d\d-\d\d).obo$')

    def __init__(self, dir_, species_name, release, ontology_release):
        assert isinstance(dir_, (str, unicode))
        assert isinstance(species_name, (str, unicode))
        assert isinstance(release, (str, unicode))
        assert isinstance(ontology_release, (str, unicode))
        
        self.dir_ = dir_.rstrip(os.sep)
        self.species_name = species_name
        self.release = release
        self.ontology_release = ontology_release

    @property
    def file_name(self):
        return u'gene_association.goa_%s.%d.gz' %(self.species_name, self.release)
    
    @property
    def path(self):
        return u'%s/%s' %(self.dir_, self.file_name)

     @property
    def file_size(self):
        return os.path.getsize(self.path)

    @property
    def ontology_file_name(self):
        return self.get_ontology_file_name(self.release)
   
    @property
    def ontology_path(self):
        return u'%s/%s' %(self.dir_, self.ontology_file_name)

    @property
    def ontology_url(self):
        """Return URL of the OBO file for this annotation data."""
        return u'http://purl.obolibrary.org/obo/go/releases/%s/go-basic.obo' \
                %(self.ontology_release)

    @staticmethod
    def get_ontology_file_name(release):
        """Get the Gene Ontology file name for a given release."""
        return u'go-basic_%s.obo' %(self.ontology_release)
    
    @staticmethod
    def get_gaf_go_release(annotation_path):
        """Extract the GO release associated with a GO annotation file.

        Parameters
        ----------
        annotation_path: str or unicode
            The path name of the GO annotation file.

        Returns
        -------
        unicode
            The release of the associated Gene Ontology file.
        """
        release = None
        with misc.smart_open_read(annotation_path, mode = 'r', encoding = 'UTF-8',
                try_gzip = True) as fh:
            for l in fh:
                if l[0] != '!':
                    break
                if l.startswith('!GO-version:'):
                    release = l.split(' ')[1]
                    break
        return release
        
    @classmethod
    def from_path(cls, path):
        assert isinstance(path, (str, unicode))
        assert os.path.isfile(path)
        
        # get the species name and release from the UniProt-GOA annotation file name
        dir_, file_name = os.path.split(path)
        m = cls.file_pat.match(file_name)
        species_name = m.group(1)
        annotation_release = int(m.group(2))
        
        # find the release of the Gene Ontology associated with the annotation file
        ontology_release = self.get_gaf_go_release(annotation_path)
        # find the file anme of the Gene Ontology file and verify its existence
        ontology_file = self.get_ontology_file_name(ontology_release)
        if not os.path.isfile(ontology_file):
            raise FileNotFoundError('Gene Ontology file "%s" not found.', ontology_file)
        
        return cls(dir_, species_name, annotation_release, ontology_release)

    @classmethod
    def scan_dir(cls, dir_):
        data = []
        dir_ = dir_.rstrip(os.sep)
        for f in os.listdir(dir_):
            path = dir_ + os.sep + f
            if os.path.isfile(path) and cls.annotation_file_pat.match(f) is not None:
                # this is a UniProt-GOA annotation file
                try:
                    data.add(cls.from_path(path))
                except FileNotFoundError:
                    # corresponding Gene Ontology file not found
                    logger.warning(('UniProt-GOA annotation file "%s" ignored, because corresponding'
                            'Gene Ontology file is missing.'), path)
        data = sorted(data, key = lambda d: [d.species_name, d.annotation_release])
        return data
