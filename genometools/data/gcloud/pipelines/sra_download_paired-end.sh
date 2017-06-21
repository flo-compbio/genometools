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

# This script donwloads SRA data, converts it to gzip'ed FASTQ(s), and uploads
# it to a bucket.

SCRIPT_DIR='{{ script_dir }}'
SRA_RUN_ACC='{{ sra_run_acc }}'
OUTPUT_DIR='{{ output_dir }}'

# store everything in /root directory
cd /root

# copy startup scripts
echo "Downloading startup scripts..."
gsutil -q cp -r "${SCRIPT_DIR}/*" ./

# install docker
chmod a+x install_docker.sh
./install_docker.sh

# download the SRA file
# echo "Downloading the .sra file..."
SRR_LOCATION="${SRA_RUN_ACC:0:6}/${SRA_RUN_ACC}/${SRA_RUN_ACC}.sra"
curl -s -O "ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/${SRR_LOCATION}"

# convert to FASTQ
echo "Converting .sra file to pair of .fastq.gz files"
docker run -v /root/:/data/ quay.io/biocontainers/sra-tools:2.8.1--0 fastq-dump -O "/data" \
        --defline-seq '@$si_$ri' --defline-qual '+' --split-files --gzip "/data/${SRA_RUN_ACC}.sra"

# upload FASTQ to bucket
echo "Uploading FASTQ files to bucket..."
gsutil -q cp "${SRA_RUN_ACC}_1.fastq.gz" "${SRA_RUN_ACC}_2.fastq.gz" "${OUTPUT_DIR}/"

# delete instance
echo "Deleting instance..."
chmod a+x delete_instance.sh
./delete_instance.sh
