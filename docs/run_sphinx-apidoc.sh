#!/bin/bash

sphinx-apidoc -e -o source/api ../genometools
rm source/api/modules.rst
