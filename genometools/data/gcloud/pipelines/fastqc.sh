#!/bin/bash

# Copyright (c) 2017 Florian Wagner
#
# This file is part of GenomeTools.
#
# GenomeTools is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License, Version 3,
# as published by the Free Software Foundation.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

# This script trims reads in a FASTQ file using trimmomatic.
# It expects a single gzip'ed FASTQ as input and outputs in the same format.

# TODO: document metadata parameters

# variables 
SCRIPT_DIR='{{ script_dir }}'
INPUT_FILE='{{ input_file }}'
OUTPUT_DIR='{{ output_dir }}'

# store everything in /root directory
cd /root

# copy startup scripts
mkdir scripts
echo "Downloading startup scripts..."
gsutil -q cp -r "${SCRIPT_DIR}/*" scripts/
chmod -R a+x scripts/*

# install docker
./scripts/install_docker.sh

# download the FASTQ file
echo "Downloading the FASTQ file from the bucket..."
gsutil -q cp "$INPUT_FILE" ./
FILE_NAME=${INPUT_FILE##*/}

# create output directory
mkdir output

# run FASTQC
docker run -u root -v "/root/:/data2/" biocontainers/fastqc fastqc \
        "/data2/${FILE_NAME}" -o /data2/output

# copy results to bucket
gsutil -q cp output/* "$OUTPUT_DIR"

{% if self_destruct %}
# delete instance
echo "Deleting instance..."
./scripts/delete_instance.sh
{% endif %}
