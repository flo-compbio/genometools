"""Package for gene set enrichment (GSE) analysis using the XL-mHG test.

The `GSEAnalysis` class performs the tests, and the results are represented
by `GSEResult` objects.
"""

from .result import GSEResult
from .analysis import GSEAnalysis
