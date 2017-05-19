#!/bin/bash
# inDrop Google Cloud pipeline

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

### PARAMETERS

# input data parameters
SCRIPT_DIR='{{ script_dir }}'
BARCODE_READ_FILE='{{ barcode_read_file }}'
MRNA_READ_FILE='{{ mrna_read_file }}'
BARCODE1_FILE='{{ barcode1_file }}'
BARCODE2_FILE='{{ barcode2_file }}'
INDEX_DIR='{{ index_dir }}'  # star index

# pipeline parameters
NUM_THREADS={{ num_threads }}
MAX_READS={{ max_reads }}

# output parameters
OUTPUT_DIR='{{ output_dir }}'

# change to root's home directory
cd /root


### PREPARATION, PART 1: Install software

# download startup scripts
echo "Downloading startup scripts from Google Storage..."
mkdir scripts
gsutil -q cp -r "${SCRIPT_DIR}/*" scripts/
chmod -R a+x scripts/*

# install compiled crcmod module
echo "Installing compiled version of crcmod module..."
./scripts/install_crcmod.sh

# install docker
echo "Installing docker..."
./scripts/install_docker.sh

# install GCC
echo "Installing GCC..."
./scripts/install_gcc.sh

# install Python 3
echo "Installing Python 3..."
./scripts/install_python.sh

# install numpy
pip3 install numpy

# install cython
pip3 install cython

# install pandas
pip3 install pandas


### PREPARATION, PART 2: Download all necessary data
mkdir data

# download input data
echo "Downloading input data from Google Storage..."
gsutil -q cp "$BARCODE1_FILE" data/barcode1_list.txt
gsutil -q cp "$BARCODE2_FILE" data/barcode2_list.txt
gsutil -q cp "$BARCODE_READ_FILE" data/barcode_reads.fastq.gz
gsutil -q cp "$MRNA_READ_FILE" data/mrna_reads.fastq.gz

# download STAR index
echo "Downloading the STAR index"
gsutil -q cp -r "${INDEX_DIR}" data/
INDEX_NAME=${INDEX_DIR##*/}


### PIPELINE STEP 1: Process the reads (identify barcodes, etc.)

# generate C file with cython
mkdir build
echo "Cythoning the barcode processing script..."
cython --embed -o build/process_barcodes.c \
        scripts/inDrop/process_barcodes.pyx

# compile and link C file to generate executable
echo "Compiling and linking the barcode processing script..."
PYTHON_INCLUDES=`python3-config --includes`
PYTHON_LIB_DIRS=`python3-config --ldflags`
PYTHON_LIBS=`python3-config --libs`
NUMPY_INCLUDE=`python3 -c "import numpy as np; print(np.get_include())"`
#gcc script.c $PYTHON_INCLUDES $PYTHON_LIB_DIRS $PYTHON_LIBS -o script
gcc build/process_barcodes.c $PYTHON_INCLUDES -I${NUMPY_INCLUDE} \
        $PYTHON_LIB_DIRS $PYTHON_LIBS \
        -o build/process_barcodes

# run the executable
echo "Running the barcode processing script..."
mkdir output_proc
./build/process_barcodes \
        -rb "data/barcode_reads.fastq.gz" -rm "data/mrna_reads.fastq.gz" \
        -b1 "data/barcode1_list.txt" -b2 "data/barcode2_list.txt" \
        -o "output_proc" \
        --max-reads $MAX_READS


### PIPELINE STEP 2: Map processed reads to the genome

# run STAR
echo "Performing alignment with STAR..."
mkdir output_star
docker run -v "/root/:/host/" quay.io/biocontainers/star:2.5.3a--0 STAR \
    --runThreadN $NUM_THREADS \
    --genomeDir "/host/data/${INDEX_NAME}" \
    --readFilesIn "/host/output_proc/processed_reads.fastq" \
    --outFileNamePrefix "/host/output_star/" \
    --outMultimapperOrder "Random" \
    --outFilterMultimapNmax 1000 \
    --outSAMmultNmax 1 \
    --outSAMtype BAM SortedByCoordinate
    # --outFilterMismatchNmax 2 \
    # --outSAMstrandField intronMotif \
    # --alignSJDBoverhangMin $SPLICE_MIN_OVERHANG \
    # --seedSearchStartLmax $SEED_START_LMAX \


### PIPELINE FINAL STEP: Upload results to Google Storage

# upload barcode counts and log file from Step 1
gsutil -q cp "output_proc/barcode_counts.tsv" "${OUTPUT_DIR}/"
gsutil -q cp "output_proc/log.txt" "${OUTPUT_DIR}/read_processing_log.txt"

# upload all STAR results from Step 2
echo "Uploading the results to Google Storage..."
gsutil -q cp "output_star/*" "${OUTPUT_DIR}/"

{% if self_destruct %}
# delete instance
echo "Deleting instance..."
./scripts/delete_instance.sh
{% endif %}
