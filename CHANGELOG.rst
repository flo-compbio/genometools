..
    Copyright (c) 2015, 2016 Florian Wagner
    
    This file is part of GenomeTools.
    
    GenomeTools is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License, Version 3,
    as published by the Free Software Foundation.
    
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.
    
    You should have received a copy of the GNU Affero General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.

Changelog
=========

Version 0.2.0 (2016-12-06)
--------------------------

- added `ontology` subpackage  

- minor bufixes and improvements

Version 0.2.1 (2016-12-16)
~~~~~~~~~~~~~~~~~~~~~~~~~~

- added function `get_protein_coding_genes` to `ensembl` subpackage

- added/expanded functionality of download functions in `misc` subpackage

- minor fixes

Version 0.2.2 (2016-12-21)
~~~~~~~~~~~~~~~~~~~~~~~~~~

- fix `six` version dependency issue with readthedocs

- minor fixes

Version 0.2.3 (2016-12-22)
~~~~~~~~~~~~~~~~~~~~~~~~~~

- added docstrings in `expression` subpackage

- minor fixes

Version 0.2.4 (2017-02-14)
~~~~~~~~~~~~~~~~~~~~~~~~~~

- rewrote `ensembl.get_protein_coding_genes()` function to use pandas for
  parsing data, store more information about the genes (including position and
  orientation), and add gene sorting feature
- added `ensembl.gdc` subpackage, which includes routines for querying the
  RESTful API of the NCI's Genomic Data Commons (https:///gdc.cancer.gov)
- minor fixes

Version 0.2.5 (2017-03-04)
~~~~~~~~~~~~~~~~~~~~~~~~~~

- fixed `ensembl_extract_protein_coding_genes.py` script
- added `ncbi.taxonomy` submodule to parse taxonomy data from NCBI FTP server
