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
OUTPUT_FILE='{{ output_file }}'
TRIM_CROP='{{ trim_crop }}'
TRIM_HEADCROP='{{ trim_headcrop }}'

# store everything in /root directory
cd /root

# copy startup scripts
echo "Downloading startup scripts..."
gsutil -q cp -r "${SCRIPT_DIR}/*" ./

# install docker
chmod a+x install_docker.sh
./install_docker.sh

# download the FASTQ file
echo "Downloading the FASTQ file from the bucket..."
gsutil -q cp "$INPUT_FILE" ./
#FILE_NAME=`basename "INPUT_FILE"
FILE_NAME=${INPUT_FILE##*/}
#TRIMMED_FILE_NAME="${INPUT_FILE##*/}_trimmed.fastq.gz"

# trim using trimmomatic
docker run -v "/root/:/data/" quay.io/biocontainers/trimmomatic:0.36--3 trimmomatic SE \
    "/data/$FILE_NAME" "/data/trimmed.fastq.gz" HEADCROP:${TRIM_HEADCROP} CROP:${TRIM_CROP}
            
echo "Uploading trimmed FASTQ file to bucket..."
gsutil -q cp "/root/trimmed.fastq.gz" "$OUTPUT_FILE"

# delete instance
echo "Deleting instance..."
chmod a+x delete_instance.sh
./delete_instance.sh
