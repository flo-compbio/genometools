# Copyright (c) 2015, 2016 Florian Wagner
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

# Package for gene set enrichment (GSE) analysis using the XL-mHG test.
# The `GSEAnalysis` class performs the tests, and the results are represented
# by `GSEResult` objects.

from .result import StaticGSEResult, RankBasedGSEResult
from .analysis import GeneSetEnrichmentAnalysis

__all__ = ['GeneSetEnrichmentAnalysis',
           'StaticGSEResult', 'RankBasedGSEResult']
