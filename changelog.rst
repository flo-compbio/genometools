..
    Copyright (c) 2015 Florian Wagner
    
    This file is part of GenomeTools.
    
    GenomeTools is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License, Version 3,
    as published by the Free Software Foundation.
    
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    
    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.


Changelog
=========

Version 1.2.2 (2016-02-26)
--------------------------
- Improved quantile normalization method (it now has a smaller memory
  footprint and handles missing data)

Version 1.2.0
-------------

- added classes for handling gene sets (see `genometools.basic`)
- added classes for handling expression data (see
  `genometools.expression`)
- reorganized package structure, created sub-packages based on data source
  (e.g., `genometools.ensembl`) or method (e.g., `genometools.rnaseq`)
- added script for finding all SRA run IDs for a given SRA experiment ID
  (``sra_find_experiment_runs.py``)
- added functions for obtaining URLs for GEO data series (see
  `genometools.geo`)
- added unicode support for all scripts

Version 1.1.0
-------------

- added documentation
- converted all tabs to four spaces
- added convenience function to configure a logger (misc.configure_logger)
  using Python's `logging` module
- added logging capabilities to the scripts extract_protein_coding_genes.py
  and extract_entrez2gene.py
- added chromosome patterns for five species to extract_protein_coding_genes.py
